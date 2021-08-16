 library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
#
#--------------------------------------------------------------
# Galaxias con line ade emision son activas: type spectral 1 ,2 ,3 (e(a),e(b) , e(c))
# El resto son galaxias pasivas: 4, 5, 6 ( k, k+a, a+k)
SS<-read.table("compact_in_gg_m3_full")      ;S<-read.table("compact_in_node_m3_full")
Gf<-read.table("compact_in_field_m3_full");  FF<-read.table("compact_in_filaments_m3_full")
HH<-read.table("../catalogos/tab_gal_gru.dat"); VG<-read.table("compact_in_voids_m3_full")


bin2=3

colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(HH)[c(1,2,3,4,5,6,7)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g")

mu_no<-S$sigv;  R_no<-S$rabs1 ; N_nodo=length(R_no)

mu_gg<-SS$sigv; R_gg<-SS$rabs1 ; N_Fof=length(R_gg);

mu_fi<-FF$sigv; R_fi<-FF$rabs1 ; N_fi=length(R_fi)

mu_cp<-Gf$sigv ; R_cp<-Gf$rabs1 ; N_cam=length(R_cp)

mu_vv<-VG$sigv ; R_vv<-VG$rabs1 ; N_vv=length(R_vv)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#seleccion de galaxias miembros
# Filamentos
fil_gr_1<-matrix(nrow=N_fi,ncol= 10)
for(i in 1:N_fi){
fil_gr_1[i,1:FF$nmi[i]]<-subset(HH$mag_g-HH$mag_r, FF$igru[i]==HH$GId)
}

fil_gr<-vector('logical',sum(FF$nmi))
h=1 ;b=FF$nmi[1]

for(i in 1:N_fi){
fil_gr[h:b] <-(fil_gr_1[i,1:FF$nmi[i]])
h=b+1 ; b=b+FF$nmi[i+1]
}
#----------------------------------------------------------------------

# Nodos
nodo_gr_1<-matrix(nrow=N_nodo,ncol= 10)
for(i in 1:N_nodo){
nodo_gr_1[i,1:S$nmi[i]]<-subset(HH$mag_g-HH$mag_r, S$igru[i]==HH$GId)
}

nodo_gr<-vector('logical',sum(S$nmi))
h=1 ;b=S$nmi[1]

for(i in 1:N_nodo){
nodo_gr[h:b] <-(nodo_gr_1[i,1:S$nmi[i]])
h=b+1 ; b=b+S$nmi[i+1]
}
#----------------------------------------------------------------------
# FOF
Fof_gr_1<-matrix(nrow=N_Fof,ncol= 10)
for(i in 1:N_Fof){
Fof_gr_1[i,1:SS$nmi[i]]<-subset(HH$mag_g-HH$mag_r, SS$igru[i]==HH$GId)
}

Fof_gr<-vector('logical',sum(SS$nmi))
h=1 ;b=SS$nmi[1]

for(i in 1:N_Fof){
Fof_gr[h:b] <-(Fof_gr_1[i,1:SS$nmi[i]])
h=b+1 ; b=b+SS$nmi[i+1]
}
#----------------------------------------------------------------------
# Campo
cam_gr_1<-matrix(nrow=N_cam,ncol= 10)
for(i in 1:N_cam){
cam_gr_1[i,1:Gf$nmi[i]]<-subset(HH$mag_g-HH$mag_r, Gf$igru[i]==HH$GId)
}

cam_gr<-vector('logical',sum(Gf$nmi))
h=1 ;b=Gf$nmi[1]

for(i in 1:N_cam){
cam_gr[h:b] <-(cam_gr_1[i,1:Gf$nmi[i]])
h=b+1 ; b=b+Gf$nmi[i+1]
}

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
# voids
vv_gr_1<-matrix(nrow=N_vv,ncol= 10)
for(i in 1:N_vv){
vv_gr_1[i,1:VG$nmi[i]]<-subset(HH$mag_g-HH$mag_r, VG$igru[i]==HH$GId)
}

vv_gr<-vector('logical',sum(VG$nmi))
h=1 ;b=VG$nmi[1]

for(i in 1:N_vv){
vv_gr[h:b] <-(vv_gr_1[i,1:VG$nmi[i]])
h=b+1 ; b=b+VG$nmi[i+1]
}

#------------------------------------------------------------------------------------------------------------

v1 <- c('Voids','Field','Filaments','FOF','Node')
v2 <- c(1:5) ; v3 <- c(1.1:4.1) ;v4 <- c(1.2:4.2)
Ngal_fil=length(fil_gr);   Ngal_cam=length(cam_gr); Ngal_fof=length(Fof_gr) ;Ngal_nodo=length(nodo_gr)  ;Ngal_vv=length(vv_gr)
field_pas<-subset(cam_gr, cam_gr >=0.7); cam_pas=length(field_pas)/length(cam_gr)
fil_pas  <-subset(fil_gr, fil_gr >=0.7); fil_pas=length(fil_pas)/length(fil_gr)
Fof_pas  <-subset(Fof_gr, Fof_gr >=0.7); Fof_pas=length(Fof_pas)/length(Fof_gr)
nodo_pas <-subset(nodo_gr,nodo_gr >=0.7); nodo_pas=length(nodo_pas)/length(nodo_gr)
vv_pas <-subset(vv_gr,vv_gr >=0.7); vv_pas=length(vv_pas)/Ngal_vv


 #Errores
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
    boot=30
    count_fil<-vector(mode="logical",boot);  count_fof<-vector(mode="logical",boot)
    count_nodo<-vector(mode="logical",boot);  count_cam<-vector(mode="logical",boot);  count_vv<-vector(mode="logical",boot)
    jj=0 ;   ff=0; ii=0 ; gg=0 ; vv=0
    r<-fil_gr ; rr<-Fof_gr ;    rrr<-nodo_gr ; rrrr<-cam_gr ;rrrrr<-vv_gr 

    theta  <- function(r){r}         ;  results<-bootstrap(r,boot,theta)
    theta2 <- function(rr){rr}      ;  results2<-bootstrap(rr,boot,theta2)
    theta3 <- function(rrr){rrr}    ;  results3<-bootstrap(rrr,boot,theta3)
    theta4 <- function(rrrr){rrrr}  ;  results4<-bootstrap(rrrr,boot,theta4)
    theta5 <- function(rrrrr){rrrrr}  ;  results5<-bootstrap(rrrrr,boot,theta5)

    for(i in 1:boot){
     for(hh in 1:Ngal_fil){
      if(results$thetastar[hh,i] >= 0.7 ){ff=ff+1}
     }
     for(hh in 1:Ngal_fof){
      if(results2$thetastar[hh,i] >= 0.7){ii=ii+1}
     }
     for(hh in 1:Ngal_nodo){
      if(results3$thetastar[hh,i] >= 0.7 ){gg=gg+1}
     }
     for(hh in 1:Ngal_cam){
      if(results4$thetastar[hh,i] >= 0.7 ){jj=jj+1}
     }

     for(hh in 1:Ngal_vv){
      if(results5$thetastar[hh,i] >= 0.7 ){vv=vv+1}
     }

     count_fil[i]  <- ff/Ngal_fil
     count_fof[i]  <- ii/Ngal_fof
     count_nodo[i] <- gg/Ngal_nodo
     count_cam[i]  <- jj/Ngal_cam
     count_vv[i]  <- vv/Ngal_vv
     ff=0; ii=0; gg=0; jj=0; vv=0
     }
    error_fil<-sd(count_fil)
    error_fof<-sd(count_fof)
    error_nodo<-sd(count_nodo)
    error_cam<-sd(count_cam)
    error_vv<-sd(count_vv)

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
naranja <- rgb(1, 101/255, 0,1)
laby=TeX('$F_{(passive)}$') 
labx=TeX('log $\\[M_*/M_\\\u0298}\\]$')


plot(NA,xaxt = "n",axes=FALSE,xlab=NA,ylab='Fraction of Galaxies',
     ylim=c(0.6,1),xlim=c(0.,5.4))
axis(side = 1,at = v2,labels = v1,tck=.02)
points(v2,c(vv_pas,cam_pas,Fof_pas,fil_pas,nodo_pas),cex=1.5,pch=c(21,20,16,17,18),col=c('black',"blue",'magenta', naranja,'red'))

error<-c(error_vv,error_cam,error_fof,error_fil,error_nodo); frac<-c(vv_pas,cam_pas,Fof_pas,fil_pas,nodo_pas)
arrows(v2,frac-error, v2,frac+error,col=c('black',"blue",'magenta', naranja,'red'),angle=90,code=3,length=0.1)


magaxis(side=1, majorn=1, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)
magaxis(side=3, majorn=1, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)
magaxis(side=4, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
laby=TeX('$\\sigma$') 
labx=TeX('$M_r$')

par(family='serif')
par(mgp=c(2.1,0.9,0))    #Posicion labels 
par(cex.lab=1.4)

legend('topleft',c('Voids','Field','Groups', "Filaments","Node"),bty="n",lty=c(0,0,0,0,0), 
       col=c('black','darkblue',"magenta",naranja,'red'),horiz=FALSE,inset=0,cex=1.2,pch=c(21,20,16,17,18))
