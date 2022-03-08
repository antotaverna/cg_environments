####
# viene del brillo_superficial.r
####
#source("propscg_vs_mabs.r")
library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
library(Hmisc)
library(fabricatr)
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

#------------------------------------------------------------------------------------------------------------



#--------
bin2=3
#--------
################################################################
result <- data.frame()
xmed<-vector() ; ymed<-vector()
ic1<-vector()  ; ic2<-vector()

func_fil <- function(arg1,arg2,arg3){
	# arg1 = Mr o deltaM
	# arg2 = props: sigma, rproy, tcr
	# arg3 = length Mr
    
    q=split_quantile(x=arg1,type=3)
    size_q1=length(arg1[q==1])
    size_q2=length(arg1[q==2])
    size_q3=length(arg1[q==3])
    
    for(i in 1:3){
    xmed[i]=mean(arg1[q==i])
    ymed[i]=mean(arg2[q==i])
    ic1[i]=boxplot(arg2[q==i], plot=F)$conf[1]
    ic2[i]=boxplot(arg2[q==i], plot=F)$conf[2]
    }

    #------------Confidence Interval
    



	#------------------BOOTSTRAP
	error_f<-vector("logical",bin2)
	boot=250; ff=0; jj=0
    count_ff<-vector(mode="logical",boot) 
    
    for(k in 1:3){
    
        r<-arg2[q==k]
        theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
    
        size_q=length(arg1[q==k])
    
        for(i in 1:boot){
	      for(hh in 1:size_q){
	          try({jj=jj+results$thetastar[hh,i]}, silent=TRUE)}
	      
             count_ff[i]<-jj/size_q
             ff=0; jj=0
          }
          error_f[k]<-sd(count_ff)
        }

        for(k in 1:bin2){
           if(error_f[k] == 0){error_f[k] = 0.02}
        }
    
        result=data.frame(xmed,ymed,error_f,ic1,ic2)
        return(result)
} #fin function
################################################################


#-------------------------------------------
#-------------------------------------------
col2alpha <- function(col, alpha) {
  col_rgb <- col2rgb(col)/255
  rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alpha)
}

#-------------------------------------------
#-------------------------------------------
func_plot <- function(arg1,arg2,cc){
	# arg1 = x
	# arg2 = y
	# arg3 = color
    return_fil=func_fil(arg1,arg2)
    xx <- return_fil[[1]]; yy <- return_fil[[2]]; erry <- return_fil[[3]]
    ic1 = return_fil[[4]]; ic2 = return_fil[[5]]


    #------------------PLOT
    points(xx,yy,type='l',col=cc,main='',lwd=2,asp=-5,xaxt='n')#,yaxt='n')
    points(xx,yy,pch=16,type='p',col=cc,main='',lwd=2,asp=-5,xaxt='n',yaxt='n',cex=1)

    #points(xx,ic1,pch=16,type='p',col=cc,main='',lwd=2,asp=-5,xaxt='n',yaxt='n',cex=1)
    #points(xx,ic2,pch=16,type='p',col=cc,main='',lwd=2,asp=-5,xaxt='n',yaxt='n',cex=1)
    polygon(c(xx,rev(xx)),c(yy+erry,rev(yy-erry)),col=col2alpha(cc,0.3),border=NA)#
    #polygon(c(xx,rev(xx)),c(ic1,rev(ic2)),col=col2alpha(cc,0.3),border=NA)#
} #fin function

#-------------------------------------------
#-------------------------------------------

#cairo_pdf("scatter_mags_quartiles.pdf")

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
limxmag=c(-22.8,-20.2)
labx=TeX('$M_{bri} - 5log_{10}h$')

laby=TeX('$\\sigma$ \\[km s$^{-1}$\\]')
plot(R_fi,sig_fi,ylim=c(120,450),xlim=limxmag,xlab='',ylab=laby,xaxt='n',type="n",log='y')#,yaxt='n')

func_plot(R_fi,sig_fi,"orange")
#----------
func_plot(R_gg,sig_gg,'magenta')
#----------
func_plot(R_no,sig_no,'red')
#----------
func_plot(R_cp,sig_cp,'darkblue')
#----------
func_plot(R_vs,sig_vs,'darkgreen')
#----------
func_plot(R_vr,sig_vr,'Deepskyblue3')

#dos filas
#legend(-22.8,490,c("Node", "Filaments",'Loose'),bty="n",lty=c(1,1,1,1,1,1), col=c('red',"orange","magenta"),horiz=FALSE,inset=0,cex=0.9,pch=c(16,16,16,16,16))
#legend(-21.65,490,c('Non-Embedded','Voids-S','Voids-R'),bty="n",lty=c(1,1,1,1,1,1), col=c('darkblue','darkgreen','Deepskyblue3'),horiz=FALSE,inset=0,cex=0.9,pch=c(16,16,16,16,16))

#una fila
legend(-22.9,220,c("Node", "Filaments",'Loose','Non-Embedded','Voids-S','Voids-R'),bty="n", col=c('red',"orange","magenta",'darkblue','darkgreen','Deepskyblue3'),horiz=FALSE,inset=0,cex=0.9,pch=c(16,16,16,16,16))

text(-20.5,150,label='A',cex=1.2)

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### Mr vs Rproy #############################
###################################################################

plot(R_fi,rp_fi,ylim=c(40,125),xlim=limxmag,xlab='',ylab=TeX('$r_p$ \\[kpc h$^{-1}$ \\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')

func_plot(R_fi,rp_fi,"orange")
#----------
func_plot(R_gg,rp_gg,'magenta')
#----------
func_plot(R_no,rp_no,'red')
#----------
func_plot(R_cp,rp_cp,'darkblue')
#----------
func_plot(R_vs,rp_vs,'darkgreen')
#----------
func_plot(R_vr,rp_vr,'Deepskyblue3')


text(-22.5,45,label='B',cex=1.2)

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)


###################################################################
######################### Mr vs tcr ###############################
###################################################################
#labx=TeX('$M_{bri}$')
plot(R_fi,tcr_fi,ylim=c(0.01,0.11),xlim=limxmag,xlab='',ylab=TeX('$H_0 \\, t_{cr}$ '),
col="orange",main='',lwd=2,asp=-5,type="n",log='y')#,xaxt='n',yaxt='n')

func_plot(R_fi,tcr_fi,"orange")
#----------
func_plot(R_gg,tcr_gg,'magenta')
#----------
func_plot(R_no,tcr_no,'red')
#----------
func_plot(R_cp,tcr_cp,'darkblue')
#----------
func_plot(R_vs,tcr_vs,'darkgreen')
#----------
func_plot(R_vr,tcr_vr,'Deepskyblue3')


text(-22.5,0.015,label='C',cex=1.2)
mtext(TeX('$M_{bri}-5log_{10}(h)$'), side = 1, cex = 0.8, line = 2.2, col = "black")

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)


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

plot(R_fi,xx_fi,ylim=c(24,26),xlim=limxmag,xlab='',ylab=TeX('$\\mu$ \\[arcsec-2\\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')

func_plot(R_fi,xx_fi,"orange")
#----------
func_plot(R_gg,xx_gg,'magenta')
#----------
func_plot(R_no,xx_no,'red')
#----------
func_plot(R_cp,xx_cp,'darkblue')
#----------
func_plot(R_vs,xx_vs,'darkgreen')
#----------
func_plot(R_vr,xx_vr,'Deepskyblue3')


legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',"orange","magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

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


plot(R_fi,xx_fi,ylim=c(0.06,0.15),xlim=limxmag,xlab='',ylab=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n",log='y')#,yaxt='n')

func_plot(R_fi,xx_fi,"orange")
#----------
func_plot(R_gg,xx_gg,'magenta')
#----------
func_plot(R_no,xx_no,'red')
#----------
func_plot(R_cp,xx_cp,'darkblue')
#----------
func_plot(R_vs,xx_vs,'darkgreen')
#----------
func_plot(R_vr,xx_vr,'Deepskyblue3')


legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',"orange","magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

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

plot(R_fi,xx_fi,ylim=c(4,17.9),xlim=limxmag,xlab='',ylab=TeX('$\\theta$ \\[arcmin\\]'),
col="orange",main='',lwd=2,asp=-5,type="n")#,xaxt='n',yaxt='n')

func_plot(R_fi,xx_fi,"orange")
#----------
func_plot(R_gg,xx_gg,'magenta')
#----------
func_plot(R_no,xx_no,'red')
#----------
func_plot(R_cp,xx_cp,'darkblue')
#----------
func_plot(R_vs,xx_vs,'darkgreen')
#----------
func_plot(R_vr,xx_vr,'Deepskyblue3')


mtext(TeX('$M_{bri}-5log_{10}(h)$'), side = 1, cex = 0.8, line = 2.2, col = "black")

legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',"orange","magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)


###################################################################
######################### DeltaM12 vs sigma #######################
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

xx_no<-delta12_no  ;yy_no<-S$sigv 
xx_gg<-delta12_gg ;yy_gg<-SS$sigv 
xx_fi<-delta12_fi ;yy_fi<-FF$sigv
xx_cp<-delta12_cp ;yy_cp<-Gf$sigv 
xx_vv<-delta12_vv ;yy_vv<-VG$sigv 
xx_vr<-delta12_vr ;yy_vr<-subset(VG$sigv,VG$tipo==0) 
xx_vs<-delta12_vs ;yy_vs<-subset(VG$sigv,VG$tipo==1) 

plot(xx_fi,yy_fi,ylim=c(100,500),xlim=c(0,2.5),xlab='',ylab=TeX('$\\sigma$ \\[km s$^{-1}$\\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')

func_plot(xx_fi,yy_fi,"orange")
#----------
func_plot(xx_gg,yy_gg,'magenta')
#----------
func_plot(xx_no,yy_no,'red')
#----------
func_plot(xx_cp,yy_cp,'darkblue')
#----------
func_plot(xx_vs,yy_vs,'darkgreen')
#----------
func_plot(xx_vr,yy_vr,'Deepskyblue3')



#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs rp ##########################
###################################################################
labx=TeX('$DeltaM_{12}$')


xx_no<-delta12_no ;yy_no<-S$rp 
xx_gg<-delta12_gg ;yy_gg<-SS$rp 
xx_fi<-delta12_fi ;yy_fi<-FF$rp
xx_cp<-delta12_cp ;yy_cp<-Gf$rp 
xx_vv<-delta12_vv ;yy_vv<-VG$rp 
xx_vr<-delta12_vr ;yy_vr<-subset(VG$rp,VG$tipo==0) 
xx_vs<-delta12_vs ;yy_vs<-subset(VG$rp,VG$tipo==1) 

return_fil=func_fil(xx_fi,yy_fi,N_fi) #count_f,prop_f,colmed_f,FL_f
plot(xx_fi,yy_fi,xlim=c(0,2.5),ylim=c(40,130),xlab='',ylab=TeX('$r_p$ \\[kpc h$^{-1}$ \\]'),
col="orange",main='',lwd=2,asp=-5,xaxt='n',type="n")#,yaxt='n')

func_plot(xx_fi,yy_fi,"orange")
#----------
func_plot(xx_gg,yy_gg,'magenta')
#----------
func_plot(xx_no,yy_no,'red')
#----------
func_plot(xx_cp,yy_cp,'darkblue')
#----------
func_plot(xx_vs,yy_vs,'darkgreen')
#----------
func_plot(xx_vr,yy_vr,'Deepskyblue3')

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

###################################################################
######################### DeltaM12 vs tcr #######################
###################################################################

labx=TeX('$\\Delta M_{12}$')
xx_no<-delta12_no  ;yy_no<-S$tcr 
xx_gg<-delta12_gg ;yy_gg<-SS$tcr
xx_fi<-delta12_fi ;yy_fi<-FF$tcr
xx_cp<-delta12_cp ;yy_cp<-Gf$tcr 
xx_vv<-delta12_vv ;yy_vv<-VG$tcr 
yy_vr<-subset(VG$tcr,VG$tipo==0)
yy_vs<-subset(VG$tcr,VG$tipo==1)


plot(xx_fi,yy_fi,ylim=c(0.01,0.15),xlim=c(0,2.5),xlab='',ylab=TeX('$H_0 \\, t_{cr}$ '),
col="orange",main='',lwd=2,asp=-5,type="n")#,xaxt='n',yaxt='n')


func_plot(xx_fi,yy_fi,"orange")
#----------
func_plot(xx_gg,yy_gg,'magenta')
#----------
func_plot(xx_no,yy_no,'red')
#----------
func_plot(xx_cp,yy_cp,'darkblue')
#----------
func_plot(xx_vs,yy_vs,'darkgreen')
#----------
func_plot(xx_vr,yy_vr,'Deepskyblue3')


mtext(expression(paste(M[2],"-",M[1])), side=1, line=2, cex=0.8)
#mtext(expression(Delta~price), side = 1, cex = 1, line = 2.2, col = "black")

#magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
#magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)

#dev.off()

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################



