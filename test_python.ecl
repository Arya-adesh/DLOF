IMPORT Python3 AS Python;
IMPORT Std.System.Thorlib;

EXPORT test_python(INTEGER neighbors=20,REAL4 contamination=0.00):= MODULE
EXPORT SKLEARN_LOF(DATASET($.datasets.anomalyLay) anomaly_ds ) := FUNCTION

test_model:=RECORD
    INTEGER SI;
    REAL SK_LOF;
    REAL SK_IS_ANOMALY;
END;

STREAMED DATASET(test_model) skl_lof(STREAMED DATASET($.datasets.anomalyLay) recs, REAL neighbors, INTEGER contamination) :=
           EMBED(Python)
    
    from sklearn.neighbors import LocalOutlierFactor
    l=[]
    for recTuple in recs:
        l.append(list(recTuple[:]));
    if(con==0.0):
        clf = LocalOutlierFactor(n_neighbors=int(K))
    else:
        clf = LocalOutlierFactor(n_neighbors=int(K),contamination=float(con))
    out= clf.fit_predict(l)
    lof_val_sk=clf.negative_outlier_factor_
    
    index=0
    
    for recTuple in recs:
        yield(index, float(lof_val_sk[index]), float(out[index]))
        index+=1
ENDEMBED;



accuracy:=skl_lof(anomaly_ds, neighbors, contamination);
return accuracy;

END;
END;
