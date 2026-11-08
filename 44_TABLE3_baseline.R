# 44_TABLE3_baseline.R
# See README.md and script_output_map.csv for inputs and outputs.
raw <- readRDS('output/phase1_reconstructed_data.rds')
raw <- raw[!is.na(raw$ID) & raw$ID != '', ]
g <- factor(raw$fenotipo_D1, levels=c('G1_free','G2_delirium','G3_coma'))
stopifnot(identical(as.integer(table(g)),c(14L,20L,11L)))
num <- function(k) {stopifnot(k %in% names(raw)); as.numeric(raw[[k]])}
out <- list()
continuous <- function(label,x){
 vals <- vapply(levels(g),function(z){v<-x[g==z];v<-v[!is.na(v)];q<-quantile(v,c(.5,.25,.75),type=6);sprintf('%.2f (%.2f-%.2f)',q[1],q[2],q[3])},'')
 out[[length(out)+1L]] <<- data.frame(Variable=label,G1=vals[1],G2=vals[2],G3=vals[3],p=kruskal.test(x~g)$p.value,row.names=NULL)
}
categorical <- function(label,x){
 vals<-vapply(levels(g),function(z)sprintf('%d/%d',sum(x[g==z],na.rm=TRUE),sum(!is.na(x[g==z]))),'')
 out[[length(out)+1L]] <<- data.frame(Variable=label,G1=vals[1],G2=vals[2],G3=vals[3],p=if(label=="Mortality") NA_real_ else fisher.test(table(g,x))$p.value,row.names=NULL)
}
cols<-c('Age (months)'='idm','PRISM IV'='PRISMIV%','VIS D1'='VISd1','VIS D5'='VISD5','Sedation duration'='sedation time','WAT-1 D1'='wat-1 d1','WAT-1 D5'='wat-1 d5','PICU stay'='TIPICU','Mechanical ventilation duration'='tvmi','Follow-up duration'='time90')
for(k in names(cols))continuous(k,num(cols[[k]]))
for(k in c('midazolam','ketamine')){
 flag<-if(k=='midazolam')'benzo' else 'cetamina';dose<-if(k=='midazolam')'dbenzo EV' else 'dceta EV'
 x<-num(dose);x[num(flag)!=1]<-NA;continuous(paste(k,'dose'),x)
}
bins<-list('Female sex'=num('sexo')==0,'Comorbidity'=num('PCCC')==1,'Underweight'=num('Nutritonal status')==1,'Sepsis'=num('sepsis')==1,'Septic shock'=num('SHOCTIP')==3,'ARDS'=num('sdra')==1,'Severe ARDS'=num('sever ards')==1,'AKI (KDIGO stage >=1)'=num('KDIGO')>=1,'Renal replacement therapy'=num('TSR')==1,'WAT-1 screening positivity'=num('Iatrogenic withdrawal syndrome (IWS)')==1,'Mortality'=num('obitooverall')==1)
for(k in names(bins))categorical(k,bins[[k]])
write.csv(do.call(rbind,out),'output/TABLE3_baseline_D1.csv',row.names=FALSE)
cat('TABLE3_baseline_D1.csv:',length(out),'rows\n')
