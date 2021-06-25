  library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------


med_ca<- c( -21.81, -20.63, -19.28754, -18.136);  cam_1<- c(0.9230, 0.7523, 0.5112360,  0.1864)
 cam_2<- c( 0.7611, 0.5239, 0.2648439,  0.0622); error_ca<-c( 0.00349, 0.009101, 0.01166, 0.00347)


med_f<- c( -22.25, -21.137, -19.91323, -18.8) ;Fil_1<- c(0.9473, 0.93103, 0.79268, 0.390)
Fil_2<- c( 0.9058, 0.74263, 0.53562, 0.244) ; error_f<-c(0.00507, 0.00960, 0.02039, 0.01561)

 med_v <-c(-21.86, -20.72, -19.506, -18.28); Voi_1<- c(0.8333, 0.7608, 0.55357,  0.3500)
Voi_2  <-c( 0.724, 0.6818, 0.5594, 0.433); med_v2<-c( -21.78, -21.224, -20.653, -20.086)
error_v<-c(0.0096, 0.0220, 0.027, 0.0159)


  med_g<- c(-22.03, -20.79, -19.54, -18.20); Gro_1<- c(0.9615, 0.9004, 0.7000, 0.3157)
  Gro_2<- c(0.9269, 0.7712, 0.5621, 0.2638); error_g<- c(0.0025, 0.0069, 0.0128, 0.0057)


rojo    <- rgb(1,0,0,1)
naranja <- rgb(1, 101/255, 0,1)
verde   <- rgb(24/255, 98/255, 24/255,1)
azul    <- rgb(0,0,1,1)
par(family="serif")
par(family="serif")
par(cex.lab=1.5)       #Tamaño labels
par(cex.axis=1.)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(4,0.4,0))  #margen labels
par(mar=c(4,4,4,4))

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------

laby=TeX('Fraction of galaxies') 
labx=TeX(' $\\M_r-5 \\, \\log{(h)}}\\]$')

  par(mfrow = c(2,2))   # dibujo 4x2
  par(mar = c(0, 0, 0,0),mgp=c(-2, 0, 0)  ,oma = c(5, 5, 5, 5))



plot(NA,xaxt = "n",axes=FALSE, ylim=c(0.1,1.05),xlim=c(-23,-17.5),xlab='',ylab=laby)
points(med_g,Gro_1,type='l',col='black',lwd=3)
points(med_g,Gro_1,type='p',col='black',lwd=3,pch=16)
points(med_g,Gro_2,type='l',col=rojo,lwd=3)
points(med_g,Gro_2,type='p',col=rojo,lwd=3,pch=17)

polygon(c(med_g,rev(med_g)),c(Gro_1+error_g,rev(Gro_1-error_g)),col=rgb(0,0,0,0.4),border=NA)
polygon(c(med_g,rev(med_g)),c(Gro_2+error_g,rev(Gro_2-error_g)),col=rgb(1, 101/255, 0,0.5),border=NA)

legend(-23,0.27,c('Gal in CG in Loose','Gal in Loose'),bty="n",lty=c(0,0),       
col=c('black',rojo),horiz=FALSE,inset=0,cex=1,pch=c(16,17))

magaxis(side=1, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE)
magaxis(side=3, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=4, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)

#---------------------------------------------------------------------------------

plot(NA,xaxt = "n",axes=FALSE, ylim=c(0.1,1.05),xlim=c(-23,-17.5),xlab='',ylab='')
points(med_f,Fil_1,type='l',col='black',lwd=2)
points(med_f,Fil_1,type='p',col='black',lwd=2,pch=16)
points(med_f,Fil_2,type='l',col=naranja,lwd=2)
points(med_f,Fil_2,type='p',col=naranja,lwd=2,pch=17)

polygon(c(med_f,rev(med_f)),c(Fil_1+error_f,rev(Fil_1-error_f)),col=rgb(0,0,0,0.4),border=NA)
polygon(c(med_f,rev(med_f)),c(Fil_2+error_f,rev(Fil_2-error_f)),col=rgb(1, 101/255, 0,0.5),border=NA)

magaxis(side=1, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=3, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=4, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)

legend(-23,0.27,c('Gal in CG in Fil','Gal in Fil'),bty="n",lty=c(0,0),    
   col=c('black',naranja),horiz=FALSE,inset=0,cex=1,pch=c(16,17))


#---------------------------------------------------------------------------------

plot(NA,xaxt = "n",axes=FALSE, ylim=c(0.05,1.05),xlim=c(-23,-17.5),xlab=labx,ylab=laby)
points(med_ca,cam_1,type='l',col='black',lwd=2)
points(med_ca,cam_1,type='p',col='black',lwd=2,pch=16)
points(med_ca,cam_2,type='l',col=azul,lwd=2)
points(med_ca,cam_2,type='p',col=azul,lwd=2,pch=17)
polygon(c(med_ca,rev(med_ca)),c(cam_1+error_ca,rev(cam_1-error_ca)),col=rgb(0,0,0,0.4),border=NA)
polygon(c(med_ca,rev(med_ca)),c(cam_2+error_ca,rev(cam_2-error_ca)),col=rgb(0,0,1,0.5),border=NA)
legend(-23,0.25,c('Gal in CG in field','Gal in field.'),bty="n",lty=c(0,0),       
col=c('black',azul),horiz=FALSE,inset=0,cex=1,pch=c(16,17))


magaxis(side=1, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE)
magaxis(side=3, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=4, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)


#---------------------------------------------------------------------------------

plot(NA,xaxt = "n",axes=FALSE, ylim=c(0.05,1.05),xlim=c(-23,-17.5),xlab=labx,ylab='')
points(med_v,Voi_1,type='l',col='black',lwd=3)
points(med_v,Voi_1,type='p',col='black',lwd=3,pch=16)
points(med_v2,Voi_2,type='l',col='gray',lwd=3)
points(med_v2,Voi_2,type='p',col='gray',lwd=3,pch=17)
polygon(c(med_v,rev(med_v)),c(Voi_1+error_v,rev(Voi_1-error_v)),col=rgb(0,0,0,0.6),border=NA)
polygon(c(med_v2,rev(med_v2)),c(Voi_2+error_v,rev(Voi_2-error_v)),col=rgb(0,0,0,0.4),border=NA)
legend(-23,0.25,c('Gal in CG in voids','Gal in voids'),bty="n",lty=c(0,0),       
col=c('black','gray'),horiz=FALSE,inset=0,cex=1,pch=c(16,17))


magaxis(side=1, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=3, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=4, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)

