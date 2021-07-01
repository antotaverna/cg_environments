#source("brillo_superficial.r")
library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
#
#--------------------------------------------------------------
# Galaxias con line ade emision son activas: type spectral 1 ,2 ,3 (e(a),e(b) , e(c))
# El resto son galaxias pasivas: 4, 5, 6 ( k, k+a, a+k)
SS<-read.table("../data/compact_in_gg_m3_full")      ;S<-read.table("../data/compact_in_node_m3_full")
Gf<-read.table("../data/compact_in_field_m3_full");  FF<-read.table("../data/compact_in_filaments_m3_full")
HH<-read.table("../data/tab_gal_gru.dat"); VG<-read.table("../data/compact_in_voids_m3_full")
naranja <- rgb(1, 101/255, 0,1)

abs_cg<-read.table("../../mabs_cg.dat")

colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(HH)[c(1,2,3,4,5,6,7)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g")
delta12_no=vector()
for(i in 1:length(S$igru)){delta12_no[i]=subset(abs_cg$V7-abs_cg$V6,abs_cg$V1==S$igru[i])}
delta12_gg=vector()
for(i in 1:length(SS$igru)){delta12_gg[i]=subset(abs_cg$V7-abs_cg$V6,abs_cg$V1==SS$igru[i])}
delta12_fi=vector()
for(i in 1:length(FF$igru)){delta12_fi[i]=subset(abs_cg$V7-abs_cg$V6,abs_cg$V1==FF$igru[i])}
delta12_cp=vector()
for(i in 1:length(Gf$igru)){delta12_cp[i]=subset(abs_cg$V7-abs_cg$V6,abs_cg$V1==Gf$igru[i])}
delta12_vv=vector()
for(i in 1:length(VG$igru)){delta12_vv[i]=subset(abs_cg$V7-abs_cg$V6,abs_cg$V1==VG$igru[i])}

#--------
bin2=3
#--------

#------------------------------------------------------------------------------------------------------------
# ------ CG in filaments-------------------------------------------------------------------------------------
result <- data.frame()

func_fil <- function(arg1,arg2,arg3){
	# arg1 = Mr
	# arg2 = props: sigma, rproy, tcr
	# arg3 = length Mr

#------------------PLOT
min_f=min(arg1); max_f=max(arg1)+0.1;   
dM=(max_f-min_f)/(bin2)

error_f<-vector("logical",bin2) ;  prop_f<-vector('logical',bin2)
error2_f<-vector("logical",bin2);  count_f<-vector('logical',bin2)
colmed_f<-vector("logical",bin2);  FL_f<-vector("logical",bin2)
a_f<-vector("logical",arg3)

count_f[1:bin2]=0;#lo agregue yo 
prop_f[1:bin2]=0
for(i in 1:arg3){
    ibin=as.integer((arg1[i]-min_f)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    count_f[k]=count_f[k]+1
    prop_f[k]=prop_f[k]+arg2[i]
    a_f[i]<- arg2[i]
   }
  
for(k in 1:bin2){
   colmed_f[k]=min_f+(dM/2.)*(2*k-1)
   FL_f[k]= (prop_f[k]/count_f[k]) #normalizado Nbin
   error2_f[k]=dM
  }

#------------------BOOTSTRAP
boot=250; ff=0; jj=0
count_ff<-vector(mode="logical",boot) 
for(k in 1:bin2){
    r<-a_f[1:count_f[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
#   theta <- function(rr){rr};results2<-bootstrap(rr,boot,theta)
    for(i in 1:boot){
      for(hh in 1:count_f[k]){
      jj=jj+results$thetastar[hh,i]
      }
      
      count_ff[i]<-jj/count_f[k] 

     ff=0; jj=0
     }
     error_f[k]<-sd(count_ff)
   }
  for(k in 1:bin2){
  if(error_f[k] == 0){error_f[k] = 0.02}
  }
#------------------

  result=data.frame(colmed_f,FL_f,error_f)
  return(result)
} #fin function
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
laby=TeX('$\\sigma$') 
labx=TeX('$M_{bri}$')

#pdf("scatter_mags.pdf")

par(family="serif")
par(cex.lab=1.5)       #Tamaño labels
par(cex.axis=0.9)      #Tamaño números en ejes.
par(lwd=2)
par(cex=2)
par(mgp=c(1.,0.2,0))   #margen labels
par(mar=c(0,3,0,1))  #c(b,l,t,r)
par(oma=c(6,1,5,1))  #c(b,l,t,r)
#par(mfrow=c(3,3))
par(mfcol=c(3,3))

R_no<-S$rabs1 ;  N_no=length(R_no)
R_gg<-SS$rabs1 ; N_gg=length(R_gg)
R_fi<-FF$rabs1 ; N_fi=length(R_fi)
R_cp<-Gf$rabs1 ; N_cp=length(R_cp)
R_vv<-VG$rabs1;  N_vv=length(R_vv)

#mu_vv<-subset(VG$tcr, VG$V15==1) ; R_vv<-subset(VG$rabs1, VG$V15==1) ; N_vv=length(R_vv)

sig_no<-S$sigv   ;rp_no<-S$rp   ;tcr_no<-S$tcr  
sig_gg<-SS$sigv  ;rp_gg<-SS$rp  ;tcr_gg<-SS$tcr
sig_fi<-FF$sigv  ;rp_fi<-FF$rp  ;tcr_fi<-FF$tcr
sig_cp<-Gf$sigv  ;rp_cp<-Gf$rp  ;tcr_cp<-Gf$tcr
sig_vv<-VG$sigv  ;rp_vv<-VG$rp  ;tcr_vv<-VG$tcr


  
###################################################################
######################### Mr vs Sigma #############################
###################################################################
 

return_fil=func_fil(R_fi,sig_fi,N_fi) #colmed_f,FL_f,error_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(100,600),xlim=c(-20,-23.),xlab='',ylab=TeX('$\\sigma$ \\[km s$^{-1}$\\]'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=16,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)
#arrows(colmed_f,FL_f-error_f, colmed_f,FL_f+error_f,col=naranja,angle=90, code=3,length=0.1)

#----------
return_gg=func_fil(R_gg,sig_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(R_no,sig_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(R_cp,sig_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(R_vv,sig_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


legend(-20.,590,c("Node", "Filaments",'Loose','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',naranja,"magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.9,pch=c(16,16,16,16,16))

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### Mr vs Rproy #############################
###################################################################

return_fil=func_fil(R_fi,rp_fi,N_fi) #colmed_f,FL_f,error_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(20,130),xlim=c(-20,-23.),xlab='',ylab=TeX('$r_p$ \\[kpc h$^{-1}$ \\]'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=16,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)
#arrows(colmed_f,FL_f-error_f, colmed_f,FL_f+error_f,col=naranja,angle=90, code=3,length=0.1)

#----------
return_gg=func_fil(R_gg,rp_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(R_no,rp_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(R_cp,rp_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(R_vv,rp_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


#legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',naranja,"magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)


###################################################################
######################### Mr vs tcr ###############################
###################################################################
#labx=TeX('$M_{bri}$')
return_fil=func_fil(R_fi,tcr_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(0,0.09),xlim=c(-20,-23.),xlab='',ylab=TeX('$H_0 \\, t_{cr}$ '),
col=naranja,main='',lwd=3,asp=-5)#,xaxt='n',yaxt='n')
points(xx,yy,pch=16,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(R_gg,tcr_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(R_no,tcr_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(R_cp,tcr_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(R_vv,tcr_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


mtext(expression(paste(M[bri])), side = 1, cex = 1, line = 2.2, col = "black")

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)


###################################################################
######################### Mr vs mu ################################
###################################################################
xx_no<-S$mu   
xx_gg<-SS$mu  
xx_fi<-FF$mu  
xx_cp<-Gf$mu  
xx_vv<-VG$mu  

return_fil=func_fil(R_fi,xx_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(23,27),xlim=c(-20,-23.),xlab='',ylab=TeX('$\\mu$ \\[arcsec-2\\]'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(R_gg,xx_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=18,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(R_no,xx_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(R_cp,xx_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=20,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(R_vv,xx_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=21,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',naranja,"magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### Mr vs dij ################################
###################################################################
xx_no<-S$dij   
xx_gg<-SS$dij  
xx_fi<-FF$dij  
xx_cp<-Gf$dij  
xx_vv<-VG$dij  


return_fil=func_fil(R_fi,xx_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(0,0.2),xlim=c(-20,-23.),xlab='',ylab=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(R_gg,xx_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=18,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(R_no,xx_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(R_cp,xx_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=20,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(R_vv,xx_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=21,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',naranja,"magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### Mr vs tita ################################
###################################################################
xx_no<-S$radio_mins*2. 
xx_gg<-SS$radio_mins*2.  
xx_fi<-FF$radio_mins*2.  
xx_cp<-Gf$radio_mins*2.  
xx_vv<-VG$radio_mins*2.  

return_fil=func_fil(R_fi,xx_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(4,15),xlim=c(-20,-23.),xlab='',ylab=TeX('$\\theta$ \\[arcmin\\]'),
col=naranja,main='',lwd=3,asp=-5)#,xaxt='n',yaxt='n')
points(xx,yy,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(R_gg,xx_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=18,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(R_no,xx_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(R_cp,xx_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=20,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(R_vv,xx_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=21,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)

mtext(expression(paste(M[bri])), side = 1, cex = 1, line = 2.2, col = "black")

legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',naranja,"magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs rp ##########################
###################################################################
bin2=3
labx=TeX('$DeltaM_{12}$')
xx_no<-delta12_no  ;yy_no<-S$rp 
xx_gg<-delta12_gg ;yy_gg<-SS$rp 
xx_fi<-delta12_fi ;yy_fi<-FF$rp
xx_cp<-delta12_cp ;yy_cp<-Gf$rp 
xx_vv<-delta12_vv ;yy_vv<-VG$rp 

return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',xlim=c(0,3),ylim=c(20,140),xlab='',ylab=TeX('$r_p$ \\[kpc h$^{-1}$ \\]'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(xx_gg,yy_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=18,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(xx_no,yy_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(xx_cp,yy_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=20,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(xx_vv,yy_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=21,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs sigma #######################
###################################################################

labx=TeX('$DeltaM_{12}$')
xx_no<-delta12_no  ;yy_no<-S$sigv 
xx_gg<-delta12_gg ;yy_gg<-SS$sigv 
xx_fi<-delta12_fi ;yy_fi<-FF$sigv
xx_cp<-delta12_cp ;yy_cp<-Gf$sigv 
xx_vv<-delta12_vv ;yy_vv<-VG$sigv 


return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(100,600),xlim=c(0,3),xlab='',ylab=TeX('$\\sigma$ \\[km s$^{-1}$\\]'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(xx_gg,yy_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=18,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(xx_no,yy_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(xx_cp,yy_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=20,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(xx_vv,yy_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=21,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs sigma #######################
###################################################################

labx=TeX('$DeltaM_{12}$')
xx_no<-delta12_no  ;yy_no<-S$tcr 
xx_gg<-delta12_gg ;yy_gg<-SS$tcr
xx_fi<-delta12_fi ;yy_fi<-FF$tcr
xx_cp<-delta12_cp ;yy_cp<-Gf$tcr 
xx_vv<-delta12_vv ;yy_vv<-VG$tcr 


return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,type='l',ylim=c(0,0.1),xlim=c(0,3),xlab='',ylab=TeX('$H_0 \\, t_{cr}$ '),
col=naranja,main='',lwd=3,asp=-5)#,xaxt='n',yaxt='n')
points(xx,yy,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_gg=func_fil(xx_gg,yy_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
points(xx,yy,type='l',col='magenta',main='',lwd=3)
points(xx,yy,pch=18,type='p',col='magenta',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

#----------
return_no=func_fil(xx_no,yy_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
points(xx,yy,pch=16,type='l',col='red',main='',lwd=3)
points(xx,yy,pch=16,type='p',col='red',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

#----------
return_cp=func_fil(xx_cp,yy_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
points(xx,yy,type='l',col='darkblue',main='',lwd=3)
points(xx,yy,pch=20,type='p',col='darkblue',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)

#----------
return_vv=func_fil(xx_vv,yy_vv,N_vv) #colmed_f,FL_f,error_f
xx <- return_vv[[1]]; yy <- return_vv[[2]]; erry <- return_vv[[3]]
points(xx,yy,type='l',col='black',main='',lwd=3)
points(xx,yy,pch=21,type='p',col='black',main='')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0,0,0,0.4),border=NA)
#arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)


mtext(expression(paste(DeltaM[12])), side = 1, cex = 1, line = 2.2, col = "black")

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

#dev.off()
