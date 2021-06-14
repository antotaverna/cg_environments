 library(bootstrap); library(latex2exp); library(astrolibR)
 library(magicaxis)
#--------------------------------------------------------------
# Galaxias con line ade emision son activas: type spectral 1 ,2 ,3 (e(a),e(b) , e(c))
# El resto son galaxias pasivas: 4, 5, 6 ( k, k+a, a+k)
SS<-read.table("compact_in_gg_m3")      ;S<-read.table("compact_in_node_m3")
Gf<-read.table("compact_in_field_m3");  FF<-read.table("compact_in_filaments_m3")

R_no<-hist(S$V14, plot=FALSE)
R_gg<-hist(SS$V11, plot=FALSE) 
R_fi<-hist(FF$V14, plot=FALSE)
R_cp<-hist(Gf$V11, plot=FALSE) 

  plot(R_no$mids,R_no$counts/sum(R_no$counts),type='s',col='red',lwd=4)
points(R_gg$mids,R_gg$counts/sum(R_gg$counts),type='s',col='magenta',lwd=4)
points(R_fi$mids,R_fi$counts/sum(R_fi$counts),type='s',col='orange',lwd=4)
points(R_cp$mids,R_cp$counts/sum(R_cp$counts),type='s',col='blue',lwd=4)
