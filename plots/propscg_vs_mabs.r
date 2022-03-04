####
# viene del brillo_superficial.r
####
#source("propscg_vs_mabs.r")
library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
library(Hmisc)
#
#
#--------------------------------------------------------------
# Galaxias con line ade emision son activas: type spectral 1 ,2 ,3 (e(a),e(b) , e(c))
# El resto son galaxias pasivas: 4, 5, 6 ( k, k+a, a+k)
SS<-read.table("../data/compact_in_gg_m3_full")      ;S<-read.table("../data/compact_in_node_m3_full")
Gf<-read.table("../data/compact_in_field_m3_full");  FF<-read.table("../data/compact_in_filaments_m3_full")
HH<-read.table("../data/tab_gal_gru.dat"); VG<-read.table("../data/compact_in_voids_m3_full")


colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')
                                                                                                                                                                          
colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')
                                                                                                                                                                          
colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <-  c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')
                                                                                                                                                                          
colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')
                                                                                                                                                                          
colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2','tipo','estado')

colnames(HH)[c(1,2,3,4,5,6,7)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g")

#--------
bin2=3
#--------
laby=TeX('$\\sigma$') 
labx=TeX('$M_{bri} - 5log_{10}h$')

#pdf("scatter_mags_eq_bins.pdf")
#pdf("scatter_mags_eq_bins_2.pdf")

par(family="serif")
par(cex.lab=1.1)       #Tamaño labels
par(cex.axis=0.8)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
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

R_vr=subset(VG$rabs1,VG$tipo==0)
R_vs=subset(VG$rabs1,VG$tipo==1)
N_vr=length(subset(VG$rabs1,VG$tipo==0))
N_vs=length(subset(VG$rabs1,VG$tipo==1))

sig_no<-S$sigv   ;rp_no<-S$rp   ;tcr_no<-S$tcr  
sig_gg<-SS$sigv  ;rp_gg<-SS$rp  ;tcr_gg<-SS$tcr
sig_fi<-FF$sigv  ;rp_fi<-FF$rp  ;tcr_fi<-FF$tcr
sig_cp<-Gf$sigv  ;rp_cp<-Gf$rp  ;tcr_cp<-Gf$tcr
sig_vv<-VG$sigv  ;rp_vv<-VG$rp  ;tcr_vv<-VG$tcr

sig_vr=subset(VG$sigv,VG$tipo==0)
sig_vs=subset(VG$sigv,VG$tipo==1)
rp_vr=subset(VG$rp,VG$tipo==0)
rp_vs=subset(VG$rp,VG$tipo==1)
tcr_vr=subset(VG$tcr,VG$tipo==0)
tcr_vs=subset(VG$tcr,VG$tipo==1)


##################################################
##################################################
##################################################
#### VOIDS S y R
##################################################
##################################################


#--------
#------------------------------------------------------------------------------------------------------------
# ------ CG in filaments-------------------------------------------------------------------------------------
R_no<-S$rabs1
R_gg<-SS$rabs1
R_fi<-FF$rabs1
R_cp<-Gf$rabs1
R_vv<-VG$rabs1

a=c(R_no,R_gg,R_fi,R_cp,R_vv)
amax=max(a)
amin=min(a)

result <- data.frame()

func_fil <- function(arg1,arg2,arg3){
	# arg1 = Mr o deltaM
	# arg2 = props: sigma, rproy, tcr
	# arg3 = length Mr
#------------------PLOT
min_f=min(arg1); max_f=max(arg1)+0.1;   
#min_f=amin; max_f=amax+0.1;   
dM=(max_f-min_f)/(bin2)
a = seq(0.0,1,0.33333)
quantile(R_fi, a, na.rm = TRUE)



library(classInt)

x <- classIntervals(dataset, 3, style = 'quantile')
x$var
x$brks

new_vec <- cut(R_fi , breaks = quantile(R_fi, a, na.rm = TRUE), include.lowest=TRUE)

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
          try({jj=jj+results$thetastar[hh,i]}, silent=TRUE)
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


#-------------------------------------------
#-------------------------------------------


func_plot <- function(xx,yy,cc){
	# arg1 = x
	# arg2 = y
	# arg3 = lim y
	# arg4 = color
#------------------PLOT
points(xx,yy,type='l',col=cc,main='',lwd=2,asp=-5,xaxt='n')#,yaxt='n')
points(xx,yy,pch=16,type='p',col=cc,main='',lwd=2,asp=-5,xaxt='n',yaxt='n',cex=1)
#polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)
} #fin function
#-------------------------------------------
#-------------------------------------------
laby=TeX('$\\sigma$') 
labx=TeX('$M_{bri} - 5log_{10}h$')

#pdf("scatter_mags_eq_bins_vRyS.pdf")
#pdf("scatter_mags_eq_bins_2_vRyS.pdf")

par(family="serif")
par(cex.lab=1.1)       #Tamaño labels
par(cex.axis=0.8)      #Tamaño números en ejes.
par(lwd=1)
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

R_vr=subset(VG$rabs1,VG$tipo==0)
R_vs=subset(VG$rabs1,VG$tipo==1)
N_vr=length(subset(VG$rabs1,VG$tipo==0))
N_vs=length(subset(VG$rabs1,VG$tipo==1))

sig_no<-S$sigv   ;rp_no<-S$rp   ;tcr_no<-S$tcr  
sig_gg<-SS$sigv  ;rp_gg<-SS$rp  ;tcr_gg<-SS$tcr
sig_fi<-FF$sigv  ;rp_fi<-FF$rp  ;tcr_fi<-FF$tcr
sig_cp<-Gf$sigv  ;rp_cp<-Gf$rp  ;tcr_cp<-Gf$tcr
sig_vv<-VG$sigv  ;rp_vv<-VG$rp  ;tcr_vv<-VG$tcr

sig_vr=subset(VG$sigv,VG$tipo==0)
sig_vs=subset(VG$sigv,VG$tipo==1)
rp_vr=subset(VG$rp,VG$tipo==0)
rp_vs=subset(VG$rp,VG$tipo==1)
tcr_vr=subset(VG$tcr,VG$tipo==0)
tcr_vs=subset(VG$tcr,VG$tipo==1)

  
###################################################################
######################### Mr vs Sigma #############################
###################################################################
 

return_fil=func_fil(R_fi,sig_fi,N_fi) #colmed_f,FL_f,error_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
laby=TeX('$\\sigma$ \\[km s$^{-1}$\\]')
plot(xx,yy,ylim=c(100,600),xlim=limxmag,xlab=c(100,600),ylab=laby,xaxt='n',type="n")#,yaxt='n')
func_plot(xx,yy,"orange")
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 0, 0,0.3),border=NA)

#----------
return_gg=func_fil(R_gg,sig_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1,0,0,0.3),border=NA)

#----------
return_no=func_fil(R_no,sig_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 101/255, 0,0.5),border=NA)

#----------
return_cp=func_fil(R_cp,sig_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0, 0, 1,0.3),border=NA)

#----------
return_vs=func_fil(R_vs,sig_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(0, 1, 0,0.3),border=NA)

#----------
return_vr=func_fil(R_vr,sig_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')
polygon(c(xx,rev(xx)),c(yy+erry*0.5,rev(yy-erry*0.5)),col=rgb(1, 1, 0,0.3),border=NA)


legend(-23.,590,c("Node", "Filaments",'Loose'),bty="n",lty=c(1,1,1,1,1,1), col=c('red',"orange","magenta"),horiz=FALSE,inset=0,cex=0.9,pch=c(16,16,16,16,16))
legend(-21.65,590,c('Non-Embedded','Voids-S','Voids-R'),bty="n",lty=c(1,1,1,1,1,1), col=c('darkblue','darkgreen','Deepskyblue3'),horiz=FALSE,inset=0,cex=0.9,pch=c(16,16,16,16,16))

text(-22.5,150,label='A',cex=1.2)

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### Mr vs Rproy #############################
###################################################################

return_fil=func_fil(R_fi,rp_fi,N_fi) #colmed_f,FL_f,error_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,ylim=c(20,130),xlim=limxmag,xlab='',ylab=TeX('$r_p$ \\[kpc h$^{-1}$ \\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(R_gg,rp_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(R_no,rp_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(R_cp,rp_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(R_vs,rp_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(R_vr,rp_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')


text(-22.5,30,label='B',cex=1.2)

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
plot(xx,yy,ylim=c(0,0.09),xlim=limxmag,xlab='',ylab=TeX('$H_0 \\, t_{cr}$ '),
col="orange",main='',lwd=2,asp=-5,type="n")#,xaxt='n',yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(R_gg,tcr_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(R_no,tcr_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(R_cp,tcr_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(R_vs,tcr_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(R_vr,tcr_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')



text(-22.5,0.01,label='C',cex=1.2)
mtext(TeX('$M_{bri}-5log_{10}(h)$'), side = 1, cex = 0.8, line = 2.2, col = "black")

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
xx_vr=subset(VG$mu,VG$tipo==0)
xx_vs=subset(VG$mu,VG$tipo==1)

return_fil=func_fil(R_fi,xx_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,ylim=c(23,27),xlim=limxmag,xlab='',ylab=TeX('$\\mu$ \\[arcsec-2\\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(R_gg,xx_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(R_no,xx_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(R_cp,xx_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(R_vs,xx_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(R_vr,xx_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')


legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',"orange","magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

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
xx_vr=subset(VG$dij,VG$tipo==0)
xx_vs=subset(VG$dij,VG$tipo==1)


return_fil=func_fil(R_fi,xx_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,ylim=c(0,0.2),xlim=limxmag,xlab='',ylab=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(R_gg,xx_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(R_no,xx_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(R_cp,xx_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(R_vs,xx_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(R_vr,xx_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')



legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',"orange","magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

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
xx_vr=subset(VG$radio_mins*2.,VG$tipo==0)
xx_vs=subset(VG$radio_mins*2.,VG$tipo==1)

return_fil=func_fil(R_fi,xx_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,ylim=c(4,15),xlim=limxmag,xlab='',ylab=TeX('$\\theta$ \\[arcmin\\]'),
col="orange",main='',lwd=2,asp=-5,type="n")#,xaxt='n',yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(R_gg,xx_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(R_no,xx_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(R_cp,xx_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(R_vs,xx_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(R_vr,xx_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')



mtext(TeX('$M_{bri}-5log_{10}(h)$'), side = 1, cex = 0.8, line = 2.2, col = "black")

legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',"orange","magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs rp ##########################
###################################################################
labx=TeX('$DeltaM_{12}$')

delta12_no= S$rabs2- S$rabs1
delta12_gg=SS$rabs2-SS$rabs1
delta12_fi=FF$rabs2-FF$rabs1
delta12_cp=Gf$rabs2-Gf$rabs1
delta12_vv=VG$rabs2-VG$rabs1
delta12_vr=subset(delta12_vv,VG$tipo==0)
delta12_vs=subset(delta12_vv,VG$tipo==1)

a=c(delta12_no,delta12_gg,delta12_fi,delta12_cp,delta12_vv)
amax=max(a)
amin=min(a)

result <- data.frame()

func_fil <- function(arg1,arg2,arg3){
	# arg1 = Mr o deltaM
	# arg2 = props: sigma, rproy, tcr
	# arg3 = length Mr
#------------------PLOT
#min_f=min(arg1); max_f=max(arg1)+0.1;   
#dM=(max_f-min_f)/(bin2)
min_f=amin; max_f=amax+0.1;   
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

#------------------

  result=data.frame(colmed_f,FL_f,error_f)
  return(result)
} #fin function



xx_no<-delta12_no  ;yy_no<-S$rp 
xx_gg<-delta12_gg ;yy_gg<-SS$rp 
xx_fi<-delta12_fi ;yy_fi<-FF$rp
xx_cp<-delta12_cp ;yy_cp<-Gf$rp 
xx_vv<-delta12_vv ;yy_vv<-VG$rp 
xx_vr<-delta12_vr ;yy_vr<-subset(VG$rp,VG$tipo==0) 
xx_vs<-delta12_vs ;yy_vs<-subset(VG$rp,VG$tipo==1) 

return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,xlim=c(0,3),ylim=c(20,140),xlab='',ylab=TeX('$r_p$ \\[kpc h$^{-1}$ \\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(xx_gg,yy_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(xx_no,yy_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(xx_cp,yy_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(xx_vs,yy_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(xx_vr,yy_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')


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
yy_vr<-subset(VG$sigv,VG$tipo==0)
yy_vs<-subset(VG$sigv,VG$tipo==1)

return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,ylim=c(100,600),xlim=c(0,3),xlab='',ylab=TeX('$\\sigma$ \\[km s$^{-1}$\\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(xx_gg,yy_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(xx_no,yy_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(xx_cp,yy_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(xx_vs,yy_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(xx_vr,yy_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')



magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs sigma #######################
###################################################################

labx=TeX('$\\Delta M_{12}$')
xx_no<-delta12_no  ;yy_no<-S$tcr 
xx_gg<-delta12_gg ;yy_gg<-SS$tcr
xx_fi<-delta12_fi ;yy_fi<-FF$tcr
xx_cp<-delta12_cp ;yy_cp<-Gf$tcr 
xx_vv<-delta12_vv ;yy_vv<-VG$tcr 
yy_vr<-subset(VG$tcr,VG$tipo==0)
yy_vs<-subset(VG$tcr,VG$tipo==1)


return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
plot(xx,yy,ylim=c(0,0.1),xlim=c(0,3),xlab='',ylab=TeX('$H_0 \\, t_{cr}$ '),
col="orange",main='',lwd=2,asp=-5,type="n")#,xaxt='n',yaxt='n')
func_plot(xx,yy,"orange")

#----------
return_gg=func_fil(xx_gg,yy_gg,N_gg) #colmed_f,FL_f,error_f
xx <- return_gg[[1]]; yy <- return_gg[[2]]; erry <- return_gg[[3]]
func_plot(xx,yy,'magenta')

#----------
return_no=func_fil(xx_no,yy_no,N_no) #colmed_f,FL_f,error_f
xx <- return_no[[1]]; yy <- return_no[[2]]; erry <- return_no[[3]]
func_plot(xx,yy,'red')

#----------
return_cp=func_fil(xx_cp,yy_cp,N_cp) #colmed_f,FL_f,error_f
xx <- return_cp[[1]]; yy <- return_cp[[2]]; erry <- return_cp[[3]]
func_plot(xx,yy,'darkblue')

#----------
return_vs=func_fil(xx_vs,yy_vs,N_vs) #colmed_f,FL_f,error_f
xx <- return_vs[[1]]; yy <- return_vs[[2]]; erry <- return_vs[[3]]
func_plot(xx,yy,'darkgreen')

#----------
return_vr=func_fil(xx_vr,yy_vr,N_vr) #colmed_f,FL_f,error_f
xx <- return_vr[[1]]; yy <- return_vr[[2]]; erry <- return_vr[[3]]
func_plot(xx,yy,'Deepskyblue3')


mtext(expression(paste(M[2],"-",M[1])), side=1, line=2, cex=0.8)
#mtext(expression(Delta~price), side = 1, cex = 1, line = 2.2, col = "black")

magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

#dev.off()

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################



