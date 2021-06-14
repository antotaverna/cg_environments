 library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
#--------------------------------------------------------------------------------
S<-read.table("compact_in_node_m3_full")
GS<-read.table("../catalogos/final_grufof_408.dat")
cg_in_fof<-read.table("../catalogos/cg_in_gg_clasifi.dat")


#--------------------------------------------------------------------------------

 ID_cg<-S$V1 ;  ID_fof<-GS$V1; N_fof
sig_cg<-S$V8 ; sig_fof<-GS$V6

flag<-cg_in_fof$V3


Id_cgin   <-intersect(cg_in_fof$V1,ID_cg )
new_ID_fof<-vector('logical',length(Id_cgin) )
jj=0
for(i in 1:length(Id_cgin)){
jj=jj+1
new_ID_fof[jj] = subset(cg_in_fof$V4,cg_in_fof$V1== Id_cgin[i])
}


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


N_cg=length(ID_cg); N_fof=length(ID_fof)
#--------------------------------------------------------------------------------
sigma1<-vector('logical', N_cg)
ll=0

 for(j in 1:47){
ll=ll+1
sigma1[ll]<-subset(sig_fof, ID_fof  == new_ID_fof[j])
}


plot(sigma1,sig_cg,pch=16)
