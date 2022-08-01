IMPORT Python3 AS Python;
IMPORT Std.System.Thorlib;
IMPORT STD;
IMPORT $;
// INTEGER C:=con* COUNT(anomaly);

EXPORT LOF(REAL contamination = 0.0,
                  UNSIGNED4 neighbors = 20,
                      STRING metric = 'minkowski',
                          STRING algorithm ='kd_tree'):= MODULE


EXPORT knn_rec4:=RECORD
    INTEGER4 SI;
    REAL4 LOF_val;
END;


EXPORT  lof(DATASET($.datasets.anomalyLay) anomaly) := FUNCTION

n_neighbors:=MIN(neighbors+1, COUNT(anomaly));
handleRec:= RECORD
  UNSIGNED handle;
END;

rec1:=RECORD  
  INTEGER4 SI;
  $.datasets.anomalyLay;
END;


// Function to initialise KD tree locally on each node
STREAMED DATASET(handleRec) fmInit(STREAMED DATASET(rec1) recs, STRING met, STRING algo) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
  
    global OBJECT
    #, OBJ_NUM
    
    import numpy 
    from sklearn.neighbors import KDTree
    from sklearn.neighbors import BallTree
        
    if 'OBJECT' not in globals():
        # This is your one-time initializer code.  It will only be executed once on each node.
        # All global initialization goes here.
       #UNPACK
        class kdCreate:
            def __init__(self,points):
                self.points=points
                l=[]
                for x in range(0, len(points)):
                    l.append(points[x][1:])
                if(algo=='kd_tree'):
                    self.tree=KDTree(l, metric=met)
                else:
                    self.tree=BallTree(l, metric=met)
                
        
        points=[] 
        for recTuple in recs:
            interList=list(recTuple[:])
            interList=list(map(float,interList))
            points.append(interList)
        OBJECT =kdCreate(points)
    return[(1,)]
ENDEMBED;

// Record to store points, its KNNs, its K distance, K distance of KNNs, LRD of its KNNs. 
knn_rec1:=RECORD
    INTEGER4 SI ; //Index of record
    INTEGER KNN_order; //KNN order number from 0 to K-1
    INTEGER4 knn; //Index of all KNNs of datapoint 
    REAL4 dis; //Distance of datapoint with its knn
    REAL4 Kdis; //Store Kdistance of datapoint of index SI
   
END;

//Function to query the created KD tree on each node

STREAMED DATASET(knn_rec1) knn(STREAMED DATASET(rec1) recs, UNSIGNED handle, INTEGER n_neighbors) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
    
    
    for recTuple in recs:
        searchItem=list(recTuple[1:])
        searchItem=list(map(float,searchItem))
        dis, ind=OBJECT.tree.query([searchItem], n_neighbors)
        
        for x in range(0, len(dis[0])):
            result=(int(recTuple[0]), int(x), int(OBJECT.points[ind[0][x]][0]) ,float(dis[0][x]), float(dis[0][n_neighbors-1]))
            yield (result)
ENDEMBED;

// 0  based numbering for SI (Index)
rec1 addSI(anomaly L, INTEGER C) := TRANSFORM
    SELF.SI:= C-1;
    SELF := L;
END;

firstDS:= PROJECT(anomaly, addSI(LEFT, COUNTER));
//Build global tree on each node hence distribute the whole dataset to each node
MyDS := DISTRIBUTE(firstDS, ALL);

handles:=fmInit(MyDS,metric, algorithm );

handle:=MIN(handles,handle);

//Query the KD tree for each point, distribute the load across all nodes
//Each node gets unique points, querying will give actual KNNs of those point

MyDS2:=DISTRIBUTE(firstDS, SI);

//Get KNN list for all points in dataset
// | SI| KNN_order|KNN| dis| Kdis
MyDS3 := knn(MyDS2, handle, n_neighbors);
// Unique records of data points with K distance 
MyDS4:=MyDS3(KNN_order=0);

 knn_rec2:=RECORD
    INTEGER4 SI;  //Index of datapoint
    INTEGER4 KNN; // index of KNN
    REAL4 reach_dis;// Store the Kdistance of Knn of point 
 END;

knn_rec2 JoinThem(MyDS3 L, MyDS4 R) := TRANSFORM
   SELF.SI:=L.SI;
   SELF.KNN:=L.KNN;
   SELF.reach_dis:= MAX(R.kdis, L.dis); //reachiblity distance  
END;

// Join to find the Kdistance of all Knns of every point
// Using Local join is much better than using global join  
withKdis:= JOIN(MyDS3,
                MyDS4,
                LEFT.KNN=RIGHT.SI and LEFT.KNN_order!=0, //Omit the point itself as Knn 
                JoinThem(LEFT, RIGHT));

// local Table to sum the reachability distance and find lrd
dis_withKdis:=DISTRIBUTE(withKdis, SI);
LRD_list := TABLE(dis_withKdis, {si, LRD := (n_neighbors-1)/SUM(GROUP, reach_dis)}, si, LOCAL);
// Local Table gave better performance than global table
// LRD_list := TABLE(withKdis, {si, LRD :=1/((SUM(GROUP, reach_dis))/ n_neighbors-1)}, si);

knn_rec3:=RECORD
    INTEGER4 SI;
    INTEGER4 KNN;
    REAL4 KNN_LRD; // LRD of all KNNs of point with index SI
 END;

knn_rec3 add_LRD(MyDS3 L, LRD_list  R) := TRANSFORM
   SELF.SI:=L.SI;
   SELF.KNN:=L.KNN;
   SELF.KNN_LRD:=R.LRD;   
END;
// Join to populate the record for the LRDs of all KNNs of every point
LRDadded:= JOIN(MyDS3,
                LRD_list,
                LEFT.KNN=RIGHT.SI and LEFT.KNN_order!=0,
                add_LRD(LEFT, RIGHT));

// Local Table gave better results than global Table
dis_LRDadded:=DISTRIBUTE(LRDadded, SI);
LOF_list := TABLE(dis_LRDadded, {si, dummy_LOF := (SUM(GROUP, KNN_LRD))/(n_neighbors-1)}, si, LOCAL);
//LOF_list := TABLE(LRDadded, {si, dummy_LOF := (SUM(GROUP, LRD2))/(K-1)}, si);
sort_LOF:=SORT(LOF_list(si>=0), si);

knn_rec4 add_LOF(sort_LOF L, LRD_list  R) := TRANSFORM
   SELF.SI:=L.SI;
   SELF.LOF_val:=-1*L.dummy_LOF/R.LRD;   
END;

//join to divide LOF value of each datapoint in record with by LRD of same point 
LOF_final:= JOIN(sort_LOF,
                LRD_list,
                LEFT.SI=RIGHT.SI,
                add_LOF(LEFT, RIGHT));

RETURN LOF_final;
END;

//////////////////////////////////////////////////////Impoved LOF /////////////////////////////////////////////////////

EXPORT  improved_lof(DATASET($.datasets.anomalyLay) anomaly) := FUNCTION

n_neighbors:=MIN(neighbors+1, COUNT(anomaly));
handleRec:= RECORD
  UNSIGNED handle;
END;

rec1:=RECORD  
  INTEGER4 SI;
  $.datasets.anomalyLay;
END;


// Function to initialise KD tree locally on each node
STREAMED DATASET(handleRec) fmInit(STREAMED DATASET(rec1) recs, STRING met, STRING algo) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
  
    global OBJECT
    #, OBJ_NUM
    
    import numpy 
    from sklearn.neighbors import KDTree
    from sklearn.neighbors import BallTree
        
    if 'OBJECT' not in globals():
        # This is your one-time initializer code.  It will only be executed once on each node.
        # All global initialization goes here.
       #UNPACK
        class kdCreate:
            def __init__(self,points):
                self.points=points
                l=[]
                for x in range(0, len(points)):
                    l.append(points[x][1:])
                if(algo=='kd_tree'):
                    self.tree=KDTree(l, metric=met)
                else:
                    self.tree=BallTree(l, metric=met)
                
        
        points=[] 
        for recTuple in recs:
            interList=list(recTuple[:])
            interList=list(map(float,interList))
            points.append(interList)
        OBJECT =kdCreate(points)
    return[(1,)]
ENDEMBED;

// Record to store points, its KNNs, its K distance, K distance of KNNs, LRD of its KNNs. 
knn_rec1:=RECORD
    INTEGER4 SI ; //Index of record
    INTEGER KNN_order; //KNN order number from 0 to K-1
    INTEGER4 knn; //Index of all KNNs of datapoint 
    REAL4 dis; //Distance of datapoint with its knn
    REAL4 Kdis; //Store Kdistance of datapoint of index SI
   
END;

knn_dup:= RECORD
    INTEGER NUM;
END;

STREAMED DATASET(knn_dup) get_duplicates(STREAMED DATASET(rec1) recs, UNSIGNED handle) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
    
    max_dup=0
    for recTuple in recs:
        searchItem=list(recTuple[1:])
        searchItem=list(map(float,searchItem))
        dup=int(OBJECT.tree.query_radius([list(searchItem)], r=0.00, count_only=True))
        if(max_dup<dup):
            max_dup=dup
    yield(max_dup)

        
ENDEMBED;

//Function to query the created KD tree on each node
STREAMED DATASET(knn_rec1) knn(STREAMED DATASET(rec1) recs, UNSIGNED handle, INTEGER K, INTEGER NORMAL) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
    
    
    for recTuple in recs:
        searchItem=list(recTuple[1:])
        searchItem=list(map(float,searchItem))
        dis, ind=OBJECT.tree.query([list(searchItem)], K)
        set_dis=set(dis[0])
        mod_dis=[]
        mod_ind=[]
        if(len(set_dis)< K):
            dis2, ind2=OBJECT.tree.query([list(searchItem)], K*NORMAL)
            set2=set(dis2[0])
            set3=list(set2)
            for x in range(0, K):
                i=dis2[0].index(set3[x])
                mod_dis.append(dis2[0][i])
                mod_ind.append(ind2[0][i])
        if(len(set_dis)==K):      
            for x in range(0, len(dis[0])):
                result=(int(recTuple[0]), int(x), int(OBJECT.points[ind[0][x]][0]),float(dis[0][x]),float(dis[0][K-1]))
                yield (result)
        else:
            for x in range(0, len(mod_dis)):
                res=(int(recTuple[0]) , int(x), int(OBJECT.points[mod_ind[x]][0]) ,float(mod_dis[x]),float(mod_dis[K-1]))
                yield(res)
ENDEMBED;

// 0  based numbering for SI (Index)
rec1 addSI(anomaly L, INTEGER C) := TRANSFORM
    SELF.SI:= C-1;
    SELF := L;
END;

firstDS:= PROJECT(anomaly, addSI(LEFT, COUNTER));
//Build global tree on each node hence distribute the whole dataset to each node
MyDS := DISTRIBUTE(firstDS, ALL);

handles:=fmInit(MyDS,metric, algorithm );

handle:=MIN(handles,handle);

//Query the KD tree for each point, distribute the load across all nodes
//Each node gets unique points, querying will give actual KNNs of those point

MyDS2:=DISTRIBUTE(firstDS, SI);
duplics:=get_duplicates(MyDS2, handle);
max_dup:=MAX(duplics, NUM);

//Get KNN list for all points in dataset
// | SI| KNN_order|KNN| dis| Kdis
MyDS3 := knn(MyDS2, handle, n_neighbors, max_dup);
// Unique records of data points with K distance 
MyDS4:=MyDS3(KNN_order=0);

 knn_rec2:=RECORD
    INTEGER4 SI;  //Index of datapoint
    INTEGER4 KNN; // index of KNN
    REAL4 reach_dis;// Store the Kdistance of Knn of point 
 END;

knn_rec2 JoinThem(MyDS3 L, MyDS4 R) := TRANSFORM
   SELF.SI:=L.SI;
   SELF.KNN:=L.KNN;
   SELF.reach_dis:= MAX(R.kdis, L.dis); //reachiblity distance  
END;

// Join to find the Kdistance of all Knns of every point
// Using Local join is much better than using global join  
withKdis:= JOIN(MyDS3,
                MyDS4,
                LEFT.KNN=RIGHT.SI and LEFT.KNN_order!=0, //Omit the point itself as Knn 
                JoinThem(LEFT, RIGHT));

// local Table to sum the reachability distance and find lrd
dis_withKdis:=DISTRIBUTE(withKdis, SI);
LRD_list := TABLE(dis_withKdis, {si, LRD := (n_neighbors-1)/SUM(GROUP, reach_dis)}, si, LOCAL);
// Local Table gave better performance than global table
// LRD_list := TABLE(withKdis, {si, LRD :=1/((SUM(GROUP, reach_dis))/ n_neighbors-1)}, si);

knn_rec3:=RECORD
    INTEGER4 SI;
    INTEGER4 KNN;
    REAL4 KNN_LRD; // LRD of all KNNs of point with index SI
 END;

knn_rec3 add_LRD(MyDS3 L, LRD_list  R) := TRANSFORM
   SELF.SI:=L.SI;
   SELF.KNN:=L.KNN;
   SELF.KNN_LRD:=R.LRD;   
END;
// Join to populate the record for the LRDs of all KNNs of every point
LRDadded:= JOIN(MyDS3,
                LRD_list,
                LEFT.KNN=RIGHT.SI and LEFT.KNN_order!=0,
                add_LRD(LEFT, RIGHT));

// Local Table gave better results than global Table
dis_LRDadded:=DISTRIBUTE(LRDadded, SI);
LOF_list := TABLE(dis_LRDadded, {si, dummy_LOF := (SUM(GROUP, KNN_LRD))/(n_neighbors-1)}, si, LOCAL);
//LOF_list := TABLE(LRDadded, {si, dummy_LOF := (SUM(GROUP, LRD2))/(K-1)}, si);
sort_LOF:=SORT(LOF_list(si>=0), si);

knn_rec4 add_LOF(sort_LOF L, LRD_list  R) := TRANSFORM
   SELF.SI:=L.SI;
   SELF.LOF_val:=-1*L.dummy_LOF/R.LRD;   
END;

//join to divide LOF value of each datapoint in record with by LRD of same point 
LOF_final:= JOIN(sort_LOF,
                LRD_list,
                LEFT.SI=RIGHT.SI,
                add_LOF(LEFT, RIGHT));

RETURN LOF_final;
END;

EXPORT without_con( DATASET(knn_rec4) LOF_final):=FUNCTION
knn_rec5:=RECORD
INTEGER SI;
REAL LOF_val;
BOOLEAN is_anomaly;
END;

// if no contamination is specified , threshold=1.5
// if LOF_val> Threshold, outlier else inlier
knn_rec5 isoutlier(LOF_final L):=TRANSFORM
   SELF.is_anomaly:=if(L.LOF_val<=-1.5, true, false);
   SELF.LOF_val:=L.LOF_val;
   SELF.SI:=L.SI;
END;
//Result with LOF value and boolean LOF value
check_anomaly:= PROJECT(LOF_final,  isoutlier(LEFT));

return check_anomaly;
END;

EXPORT with_con( DATASET(knn_rec4) LOF_final,DATASET($.datasets.anomalyLay) anomaly ):=FUNCTION
//Obtain threshold value based on contamination
INTEGER C:=contamination* COUNT(anomaly);
lofTOPN:=TOPN(LOF_final,C+1,-LOF_val);
REAL4 t:=lofTOPN[C+1].LOF_val;

knn_rec5:=RECORD
INTEGER SI;
REAL LOF_val;
BOOLEAN is_anomaly;
END;

// if LOF_val> Threshold, outlier else inlier
knn_rec5 isoutlier(LOF_final L):=TRANSFORM
   SELF.is_anomaly:=if(L.LOF_val<=t, true, false);
   SELF.LOF_val:=L.LOF_val;
   SELF.SI:=L.SI;
END;
//Result with LOF value and boolean LOF value
check_anomaly:= PROJECT(LOF_final,  isoutlier(LEFT));
return check_anomaly;
END;

EXPORT  fit_predict(DATASET(knn_rec4) LOF_final, DATASET($.datasets.anomalyLay) anomaly) := FUNCTION
out:=if(contamination=0.0, without_con(LOF_final), with_con(LOF_final, anomaly));
return out;
END;
END;









