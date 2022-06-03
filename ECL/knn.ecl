IMPORT Python3 AS Python;
anomaly:= $.File_dlof.File;
anomalyLay:=$.File_dlof.Layout;

handleRec := RECORD
  UNSIGNED handle;
END;
dummy_rec:=RECORD  
  INTEGER4 SI;
  anomalyLay;
END;
STREAMED DATASET(handleRec) fmInit(STREAMED DATASET(dummy_rec) recs) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
  
    global OBJECT
    #, OBJ_NUM
    
    import numpy 
    from sklearn.neighbors import KDTree

    if 'OBJECT' not in globals():
        # This is your one-time initializer code.  It will only be executed once on each node.
        # All global initialization goes here.
       #UNPACK
        class kdCreate:
            def __init__(self,points):
                self.kdis=[]
                self.tree=KDTree(points)
            def storeKdis(self, ind, dis):
                self.kdis.append((ind,dis))
                
        
        points=[] 
        for recTuple in recs:
            interList=list(recTuple[1:])
            interList= list(map(float,interList))
            points.append(interList)

    # Now instantiate the object that we want to use repeatedly
        OBJECT =kdCreate(points)
        

    # We return a single dummy record with the object handle inside.
    return[(1,)]
ENDEMBED;

// Here's a routine that uses the shared object from Init.
// Notice that it must receive handle even though it's not used.
// Otherwise, we can't guarantee that fmInit will be called first.
knn_rec:=RECORD
    INTEGER4 SI ;
    INTEGER4 knn;
    REAL4 dis;
    REAL4 reach:=0;
    REAL4 LRD:=0;
END;

STREAMED DATASET(knn_rec) knn(STREAMED DATASET(dummy_rec) recs, UNSIGNED handle, INTEGER K) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
    
    
    for recTuple in recs:
        searchItem=list(recTuple[1:])
        searchItem=list(map(float,searchItem))
        dis, ind=OBJECT.tree.query([list(searchItem)], K)
        
        for x in range(0, len(dis[0])):
            result=(int(recTuple[0]),int(ind[0][x]),float(dis[0][x]),float(0),float(0))
            yield (result)
        OBJECT.kdis.append((int(recTuple[0])))
ENDEMBED;

knn_rec2:= RECORD
    INTEGER SI;
    //REAL4 KDIS;
END;

STREAMED DATASET(knn_rec2) kdistance(STREAMED DATASET(knn_rec) recs,UNSIGNED handle ) :=
           EMBED(Python: globalscope('facScope'), persist('query'), activity)
    
    
    for x in range(0,len(OBJECT.kdis)):
       yield(OBJECT.kdis[x][0])
        

ENDEMBED;

dummy_rec addSI(anomaly L, INTEGER C) := TRANSFORM
    SELF.SI:= C-1;
    SELF := L;
END;

firstDS:= PROJECT(anomaly, addSI(LEFT, COUNTER));
MyDS := DISTRIBUTE(firstDS, ALL);
OUTPUT(MyDS, NAMED('InputDataset'));

handles:=fmInit(MyDS);
OUTPUT(handles, NAMED('handles'));
handle:=MIN(handles,handle);
OUTPUT(handle, NAMED('handle'));

MyDS2:=DISTRIBUTE(firstDS);////////////////////////////////////////IMPORTANT 
OUTPUT(MyDS2, NAMED('MyDS2'));
INTEGER K:=5;
MyDS3 := knn(firstDS, handle, K);
OUTPUT(MyDS3, NAMED('MyDS3'));

knn_rec3:=RECORD
    knn_rec.SI;
    knn_rec.knn;
    knn_rec.dis;
    knn_rec.reach;
    knn_rec.LRD;
END;

knn_rec JoinThem(MyDS3 L, MyDS3 R) := TRANSFORM
    SELF.SI:=L.SI;
    SELF.knn:=L.knn;
    SELF.dis:=L.dis;
    SELF.reach:=R.dis;
    SELF.LRD:=L.LRD;
   
    
END;

//d1:=dataset(MyDS3, knn_rec);

RT_folk := JOIN(MyDS3,
                MyDS3,
                LEFT.knn=RIGHT.SI,
                JoinThem(LEFT, RIGHT));
OUTPUT(RT_folk, NAMED('JOINED'));
