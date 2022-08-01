IMPORT $;
LOF_OUT:=$.LOF(neighbors:=10, contamination:=0.5).lof($.datasets.anomaly);
FIT_OUT:=$.LOF(neighbors:=10, contamination:=0.5).fit_predict(LOF_OUT,$.datasets.anomaly);
ANOMALIES:=FIT_OUT(is_anomaly=True);
OUTPUT(LOF_OUT, NAMED('LOF'));
OUTPUT(FIT_OUT, NAMED('Fit_Predict'));
OUTPUT(ANOMALIES, NAMED('Anomalies'));

