 library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
library(cosmoFns)
#
#--------------------------------------------------------------
# Galaxias con line ade emision son activas: type spectral 1 ,2 ,3 (e(a),e(b) , e(c))
# El resto son galaxias pasivas: 4, 5, 6 ( k, k+a, a+k)
  FF<-read.table("../data/compact_in_field_m3_full")
HH<-read.table("../data/tab_gal_gru.dat")

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')


colnames(HH)[c(1,2,3,4,5,6,7)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g")

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#k-corrections ([Chilingarian2010])
coef_gr<-matrix(nrow=6, ncol=4)
coef_gr[1,1:4]<- c(  0,0,0,0 )
coef_gr[2,1:4]<- c( 1.83285,-2.71446,4.97336,-3.66864 )
coef_gr[3,1:4]<- c( -19.7595,10.5033,18.8196,6.07785)
coef_gr[4,1:4]<- c( 33.6059,-120.713,-49.299,0)
coef_gr[5,1:4]<- c( 144.371,216.453,0,0 )
coef_gr[6,1:4]<- c(-295.39,0,0,0)


coefg_gr<-matrix(nrow=8, ncol=4)
coefg_gr[1,1:4]<- c(  0,0,0,0 )
coefg_gr[2,1:4]<- c(-2.45204,4.10188,10.5258,-13.5889 )
coefg_gr[3,1:4]<- c( 56.7969,-140.913,144.572,57.2155 )
coefg_gr[4,1:4]<- c(-466.949,222.789,-917.46,-78.0591 )
coefg_gr[5,1:4]<- c(2906.77,1500.8,1689.97,30.889  )
coefg_gr[6,1:4]<- c( -10453.7,-4419.56,-1011.01,0  )
coefg_gr[7,1:4]<- c( 17568,3236.68,0,0 )
coefg_gr[8,1:4]<- c(  -10820.7,0,0,0  )

              
K_r<-vector('logical',length(HH$Redshift)) ; K_g<-vector('logical',length(HH$Redshift))
K_r[1:length(HH$Redshift)]=0 ; K_g[1:length(HH$Redshift)]=0
color_gr<-HH$mag_g-HH$mag_r

for(l in 1:length(HH$Redshift)){
 for (i in 1:6){
  for (j in 1:4){
  K_r[l]= K_r[l]+ coef_gr[i,j]*((HH$Redshift[l])**(i-1)) * (color_gr[l])**(j-1)
  }
 } 
}

for(l in 1:length(HH$Redshift)){
 for (i in 1:8){
  for (j in 1:4){
  K_g[l]= K_g[l]+ coefg_gr[i,j]*((HH$Redshift[l])**(i-1)) * (color_gr[l])**(j-1)
  }
 } 
}

#Magnitudes absolutas
h_l=0.673
Mag_r<-vector('logical',length(HH$Redshift))
for(l in 1:length(HH$Redshift)){

 dL<-D.L(HH$Redshift[l], omega.m =0.315, omega.lambda =0.685, H.0 = 67.3)

 Mag_r[l]= HH$mag_r[l]-25.-5.*log10(dL)-K_r[l]
Mag_r[l]=Mag_r[l]-5*log10(h_l)
}

Mag_g<-vector('logical',length(HH$Redshift))
for(l in 1:length(HH$Redshift)){
 dL<-D.L(HH$Redshift[l], omega.m =0.315, omega.lambda =0.685, H.0 = 67.3)
 Mag_g[l]= HH$mag_g[l]-25.-5.*log10(dL)-K_g[l]
Mag_g[l]=Mag_g[l]-5*log10(h_l)
}

#------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
#seleccion de galaxias miembros de grupos compactos en
# grupos  loose
 N_fi=length(FF$rabs1)
fil_r<-matrix(nrow=N_fi,ncol= 10)
fil_g<-matrix(nrow=N_fi,ncol= 10)
for(i in 1:N_fi){
fil_r[i,1:FF$nmi[i]]<-subset(Mag_r, FF$igru[i]==HH$GId)
fil_g[i,1:FF$nmi[i]]<-subset(Mag_g, FF$igru[i]==HH$GId)
}

Mag_Cr<-vector('logical',sum(FF$nmi)); Mag_Cg<-vector('logical',sum(FF$nmi))

h=1 ;b=FF$nmi[1]

for(i in 1:N_fi){
Mag_Cr[h:b] <-(fil_r[i,1:FF$nmi[i]])
Mag_Cg[h:b] <-(fil_g[i,1:FF$nmi[i]])
h=b+1 ; b=b+FF$nmi[i+1]

}


#---------------------------------------------------------------------------------------------------------------------
# Lectura de galaxias en filamentos

  G_FF<-read.table("../data/field_galaxy_sample.dat")

colnames(G_FF)[c(1,2,3,4,5,6)] <- c("igru",'ra','dec','z','r','g')

r_apr<-subset(G_FF$r, G_FF$r > 1. )
g_apr<-subset(G_FF$g, G_FF$r > 1. )
z_apr<-subset(G_FF$z, G_FF$r > 1. )

K_r2<-vector('logical',length(z_apr)) ; K_r2[1:length(z_apr)]=0 
K_g2<-vector('logical',length(z_apr)) ; K_g2[1:length(z_apr)]=0


color_gr2<-g_apr-r_apr

for(l in 1:length(r_apr)){
 for (i in 1:6){
  for (j in 1:4){
  K_r2[l]= K_r2[l]+ coef_gr[i,j]*((z_apr[l])**(i-1)) * (color_gr2[l])**(j-1)
  }
 } 
}

for(l in 1:length(r_apr)){
 for (i in 1:8){
  for (j in 1:4){
  K_g2[l]= K_g2[l]+ coefg_gr[i,j]*((z_apr[l])**(i-1)) * (color_gr2[l])**(j-1)
  }
 } 
}


#Magnitudes absolutas
h_l=0.673
M_r<-vector('logical',length(z_apr)) ; M_r[1:length(z_apr)] =0
M_g<-vector('logical',length(z_apr)) ; M_g[1:length(z_apr)] =0


for(l in 1:length(z_apr)){
 dL<-D.L(z_apr[l], omega.m =0.315, omega.lambda =0.67, H.0 = 67.3)

 M_r[l]= r_apr[l]-25.-5.*log10(dL)-K_r2[l]
 M_r[l]=M_r[l]-5*log10(h_l)
 M_g[l]= g_apr[l]-25.-5.*log10(dL)-K_g2[l]
 M_g[l]=M_g[l]-5*log10(h_l)
}




#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
min_f=min(Mag_Cr);  max_f=max(Mag_Cr);bin2=4
dM=(max_f-min_f)/(bin2)

 count_f<-vector('logical',bin2); total_f<-vector('logical',bin2); med_f<-vector('logical',bin2)
 count_f[1:bin2]=0; total_f[1:bin2]=0
 N_fi=length(Mag_Cr)
 Mag_Cg_Mag_Cr<-matrix(nrow=N_fi,ncol=bin2)

  for(i in 1:N_fi){
   ibin=as.integer((Mag_Cr[i]-min_f)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    if(Mag_Cg[i]-Mag_Cr[i] >=0.7){count_f[k]=count_f[k]+1}
    total_f[k]=total_f[k]+1
    med_f[k]= med_f[k]+Mag_Cr[i]
    Mag_Cg_Mag_Cr[total_f[k],k]<-Mag_Cg[i]-Mag_Cr[i]
   }

FL_fil<-vector("logical",bin2); med_1<-vector("logical",bin2);
  for(i in 1:bin2){
 FL_fil[i]=count_f[i]/total_f[i]
 med_1[i]=med_f[i]/total_f[i]
}

#---------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------

    boot=300
    count_fil<-vector(mode="logical",boot);  count_fof<-vector(mode="logical",boot)
    count_nodo<-vector(mode="logical",boot);  count_cam<-vector(mode="logical",boot);  count_vv<-vector(mode="logical",boot)
    error_fil<-vector(mode='logical', bin2)
   
   for(k in 1:bin2){
    r<- Mag_Cg_Mag_Cr[1:total_f[k],k]
     theta  <- function(r){r} ;  results<-bootstrap(r,boot,theta)
    ff=0; ii=0 ;
     for(i in 1:boot){
      for(hh in 1:total_f[k]){
       if(results$thetastar[hh,i] >= 0.7 ){ff=ff+1}
      }
      count_fil[i]  <- ff/N_fi
      ff=0
     }
    error_fil[k]<-sd(count_fil)
   }



#---------------------------------------------------------------------------------------------------------------------
M_gtol<-subset(M_g, M_g > -30 &  M_r >  min(Mag_Cr) &  M_r < max(Mag_Cr) &  M_g < max(Mag_Cg) &  M_g > min(Mag_Cg))
M_rtol<-subset(M_r, M_g > -30 &  M_r >  min(Mag_Cr) &  M_r < max(Mag_Cr) &  M_g < max(Mag_Cg) &  M_g > min(Mag_Cg))
abs_color<- (M_gtol -M_rtol)

min_f=min(M_rtol);  max_f=max(M_rtol);bin2=4
dM=(max_f-min_f)/(bin2)

 count_2<-vector('logical',bin2); total_f<-vector('logical',bin2); med_f<-vector('logical',bin2)
 count_2[1:bin2]=0; total_f[1:bin2]=0
 N_fi2=length(G_FF$z)
  for(i in 1:length(M_rtol)){
   ibin=as.integer((M_rtol[i]-min_f)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    if(abs_color[i] >=0.7){count_2[k]=count_2[k]+1}
    total_f[k]=total_f[k]+1
    med_f[k]= med_f[k]+M_rtol[i]
   }


FL_2<-vector("logical",bin2); med_2<-vector("logical",bin2)

  for(i in 1:bin2){
 FL_2[i]=count_2[i]/total_f[i]
 med_2[i]=med_f[i]/total_f[i]
}




#-------------------------------------------------------------------------------------------------------
par(family="serif")
par(cex.lab=0.7)       #Tamaño labels
par(cex.axis=0.6)      #Tamaño números en ejes.
par(lwd=2)
par(cex=2)
par(mgp=c(1.,0.1,0))      	        #margen labels
par(mar=c(4,4,4,4))
##-------------------------------------------------------------
rojo    <- rgb(1,0,0,1)
naranja <- rgb(1, 101/255, 0,1)
verde   <- rgb(24/255, 98/255, 24/255,1)
azul    <- rgb(0,0,1,1)
#--------------------------------------------------------------


laby=TeX('Fraccion of galaxies') 
labx=TeX(' $\\M_r-5 \\, \\log{(h)}}\\]$')


plot(NA,xaxt = "n",axes=FALSE, ylim=c(0.,1.2),xlim=c(-23,-17.5),xlab=labx,ylab=laby)
points(med_1,FL_fil,type='l',col='black',lwd=5)
points(med_1,FL_fil,type='p',col='black',lwd=1,pch=16)
points(med_2,FL_2,type='l',col=azul,lwd=5)
points(med_2,FL_2,type='p',col=azul,lwd=1,pch=16)

polygon(c(med_1,rev(med_1)),c(FL_fil+error_fil,rev(FL_fil-error_fil)),col=rgb(0,0,0,0.4),border=NA)
polygon(c(med_2,rev(med_2)),c(FL_2+error_fil*0.6,rev(FL_2-error_fil*0.6)),col=rgb(0,0,1,0.5),border=NA)


legend(-23,0.2,c('Gal in CG in field','Gal in field.'),bty="n",lty=c(0,0),       
col=c('black',rojo),horiz=FALSE,inset=0,cex=0.6,pch=c(16,16,16,16,16))

magaxis(side=1, majorn=5, minorn=5, tcl=0.3, ratio=0.5, labels=TRUE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.3, ratio=0.5, labels=TRUE)
magaxis(side=3, majorn=5, minorn=5, tcl=0.3, ratio=0.5, labels=FALSE)
magaxis(side=4, majorn=5, minorn=5, tcl=0.3, ratio=0.5, labels=FALSE)






