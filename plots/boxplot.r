#source("boxplot.r") 
library(bootstrap); library(latex2exp); library(astrolibR)
 #library(magicaxis)
#--------------------------------------------------------------

SS<-read.table("../data/compact_in_gg_m3_full")
S<-read.table("../data/compact_in_node_m3_full")
Gf<-read.table("../data/compact_in_field_m3_full")
FF<-read.table("../data/compact_in_filaments_m3_full")
VG<-read.table("../data/compact_in_voids_m3_full")
HH<-read.table("../data/tab_gal_gru.dat")

#---------COLNAMES ---------
colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp','rabs2','tipo','estado')

colnames(HH)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g", "mag_g", "galid", "r50", "r90", "rabs", "gabs", "uabd", "passive", "early")

#R_no<-hist(S$V14, plot=FALSE)
#R_gg<-hist(SS$V11, plot=FALSE) 
#R_fi<-hist(FF$V14, plot=FALSE)
#R_cp<-hist(Gf$V11, plot=FALSE) 

#  plot(R_no$mids,R_no$counts/sum(R_no$counts),type='s',col='red',lwd=4)
#points(R_gg$mids,R_gg$counts/sum(R_gg$counts),type='s',col='magenta',lwd=4)
#points(R_fi$mids,R_fi$counts/sum(R_fi$counts),type='s',col='orange',lwd=4)
#points(R_cp$mids,R_cp$counts/sum(R_cp$counts),type='s',col='blue',lwd=4)


#------------------------------------------------------------------
#------------------------------------------------------------------
#seleccion de galaxias miembros
#------------------------------------------------------------------

N_nodo=length(S$igru)
N_Fof=length(SS$igru)
N_fi=length(FF$igru)
N_cam=length(Gf$igru)
N_vv=length(VG$igru)
N_vr=length(subset(VG$igru,VG$tipo==0))
N_vs=length(subset(VG$igru,VG$tipo==1))
N_ved=length(subset(VG$igru,VG$estado==0))
N_vin=length(subset(VG$igru,VG$estado==-1))

#------------------------------------------------------------------
# Filamentos
fil_gr_1<-matrix(nrow=N_fi,ncol= 10)
fil_frac_p1<-matrix(nrow=N_fi,ncol= 10)
fil_frac_e1<-matrix(nrow=N_fi,ncol= 10)

for(i in 1:N_fi){
    fil_gr_1[i,1:FF$nmi[i]]<-subset(HH$mag_r, FF$igru[i]==HH$GId)
    fil_frac_p1[i,1:FF$nmi[i]]<-subset(HH$passive, FF$igru[i]==HH$GId)
    fil_frac_e1[i,1:FF$nmi[i]]<-subset(HH$early, FF$igru[i]==HH$GId)
}

fil_gr<-vector('logical',sum(FF$nmi)); fil_dif_2bri<-vector('logical',N_fi)
fil_fp<-vector('logical',N_fi); fil_fp_gr<-vector('logical',N_fi)
fil_fe<-vector('logical',N_fi); fil_fe_gr<-vector('logical',N_fi)

h=1 ;b=FF$nmi[1]

for(i in 1:N_fi){
    fil_gr[h:b] <-(fil_gr_1[i,1:FF$nmi[i]]) #concatena miembros
    fil_fp[i] <-sum(fil_frac_p1[i,1:FF$nmi[i]])
    fil_fe[i] <-sum(fil_frac_e1[i,1:FF$nmi[i]])
    h=b+1 ; b=b+FF$nmi[i+1]

}

h=1 ;b=FF$nmi[1]

for(i in 1:N_fi){
    fil_dif_2bri[i]<-fil_gr[h+1]-fil_gr[h]
    fil_fp_gr[i]<-fil_fp[i]/FF$nmi[i]
    fil_fe_gr[i]<-fil_fe[i]/FF$nmi[i]
    h=b+1 ; b=b+FF$nmi[i+1]
}



#----------------------------------------------------------------------

# Nodos
nodo_gr_1<-matrix(nrow=N_nodo,ncol= 10)
nodo_frac_p1<-matrix(nrow=N_nodo,ncol= 10)
nodo_frac_e1<-matrix(nrow=N_nodo,ncol= 10)

for(i in 1:N_nodo){
    nodo_gr_1[i,1:S$nmi[i]]<-subset(HH$mag_r, S$igru[i]==HH$GId)
    nodo_frac_p1[i,1:S$nmi[i]]<-subset(HH$passive, S$igru[i]==HH$GId)
    nodo_frac_e1[i,1:S$nmi[i]]<-subset(HH$early, S$igru[i]==HH$GId)
}

nodo_gr<-vector('logical',sum(S$nmi)); nodo_dif_2bri<-vector('logical',N_nodo)
nodo_fp<-vector('logical',N_nodo); nodo_fp_gr<-vector('logical',N_nodo)
nodo_fe<-vector('logical',N_nodo); nodo_fe_gr<-vector('logical',N_nodo)

h=1 ;b=S$nmi[1]
for(i in 1:N_nodo){
    nodo_gr[h:b] <-(nodo_gr_1[i,1:S$nmi[i]])
    nodo_fp[i] <-sum(nodo_frac_p1[i,1:S$nmi[i]])
    nodo_fe[i] <-sum(nodo_frac_e1[i,1:S$nmi[i]])
    h=b+1 ; b=b+S$nmi[i+1]
}

h=1 ;b=S$nmi[1]
for(i in 1:N_nodo){
    nodo_dif_2bri[i]<-nodo_gr[h+1]-nodo_gr[h]
    nodo_fp_gr[i]<-nodo_fp[i]/S$nmi[i]
    nodo_fe_gr[i]<-nodo_fe[i]/S$nmi[i]
    h=b+1 ; b=b+S$nmi[i+1]
}

#----------------------------------------------------------------------
# FOF
Fof_gr_1<-matrix(nrow=N_Fof,ncol= 10)
Fof_frac_p1<-matrix(nrow=N_Fof,ncol= 10)
Fof_frac_e1<-matrix(nrow=N_Fof,ncol= 10)

for(i in 1:N_Fof){
    Fof_gr_1[i,1:SS$nmi[i]]<-subset(HH$mag_r, SS$igru[i]==HH$GId)
    Fof_frac_p1[i,1:SS$nmi[i]]<-subset(HH$passive, SS$igru[i]==HH$GId)
    Fof_frac_e1[i,1:SS$nmi[i]]<-subset(HH$early, SS$igru[i]==HH$GId)
}

Fof_gr<-vector('logical',sum(SS$nmi)); Fof_dif_2bri<-vector('logical',N_Fof)
Fof_fp<-vector('logical',N_Fof); Fof_fp_gr<-vector('logical',N_Fof)
Fof_fe<-vector('logical',N_Fof); Fof_fe_gr<-vector('logical',N_Fof)

h=1 ;b=SS$nmi[1]
for(i in 1:N_Fof){
    Fof_gr[h:b] <-(Fof_gr_1[i,1:SS$nmi[i]])
    Fof_fp[i] <-sum(Fof_frac_p1[i,1:SS$nmi[i]])
    Fof_fe[i] <-sum(Fof_frac_e1[i,1:SS$nmi[i]])
    h=b+1 ; b=b+SS$nmi[i+1]
}

h=1 ;b=SS$nmi[1]
for(i in 1:N_Fof){
    Fof_dif_2bri[i]<-Fof_gr[h+1]-Fof_gr[h]
    Fof_fp_gr[i]<-Fof_fp[i]/SS$nmi[i]
    Fof_fe_gr[i]<-Fof_fe[i]/SS$nmi[i]
    h=b+1 ; b=b+SS$nmi[i+1]
}


#----------------------------------------------------------------------
# Campo
cam_gr_1<-matrix(nrow=N_cam,ncol= 10)
cam_frac_p1<-matrix(nrow=N_cam,ncol= 10)
cam_frac_e1<-matrix(nrow=N_cam,ncol= 10)
for(i in 1:N_cam){
    cam_gr_1[i,1:Gf$nmi[i]]<-subset(HH$mag_r, Gf$igru[i]==HH$GId)
    cam_frac_p1[i,1:Gf$nmi[i]]<-subset(HH$passive, Gf$igru[i]==HH$GId)
    cam_frac_e1[i,1:Gf$nmi[i]]<-subset(HH$early, Gf$igru[i]==HH$GId)
}

cam_gr<-vector('logical',sum(Gf$nmi)); cam_dif_2bri<-vector('logical',N_cam)
cam_fp<-vector('logical',N_cam); cam_fp_gr<-vector('logical',N_cam)
cam_fe<-vector('logical',N_cam); cam_fe_gr<-vector('logical',N_cam)

h=1 ;b=Gf$nmi[1]
for(i in 1:N_cam){
    cam_gr[h:b] <-(cam_gr_1[i,1:Gf$nmi[i]])
    cam_fp[i] <-sum(cam_frac_p1[i,1:Gf$nmi[i]])
    cam_fe[i] <-sum(cam_frac_e1[i,1:Gf$nmi[i]])
    h=b+1 ; b=b+Gf$nmi[i+1]
}

h=1 ;b=Gf$nmi[1]
for(i in 1:N_cam){
    cam_dif_2bri[i]<-cam_gr[h+1]-cam_gr[h]
    cam_fp_gr[i]<-cam_fp[i]/Gf$nmi[i]
    cam_fe_gr[i]<-cam_fe[i]/Gf$nmi[i]
    h=b+1 ; b=b+Gf$nmi[i+1]
}


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
# voids
vv_gr_1<-matrix(nrow=N_vv,ncol= 10)
vv_frac_p1<-matrix(nrow=N_vv,ncol= 10)
vv_frac_e1<-matrix(nrow=N_vv,ncol= 10)

for(i in 1:N_vv){
    vv_gr_1[i,1:VG$nmi[i]]<-subset(HH$mag_r, VG$igru[i]==HH$GId)
    vv_frac_p1[i,1:VG$nmi[i]]<-subset(HH$passive, VG$igru[i]==HH$GId)
    vv_frac_e1[i,1:VG$nmi[i]]<-subset(HH$early, VG$igru[i]==HH$GId)
}

vv_gr<-vector('logical',sum(VG$nmi)); vv_dif_2bri<-vector('logical',N_vv)
vv_fp<-vector('logical',N_vv); vv_fp_gr<-vector('logical',N_vv)
vv_fe<-vector('logical',N_vv); vv_fe_gr<-vector('logical',N_vv)

h=1 ;b=VG$nmi[1]
for(i in 1:N_vv){
    vv_gr[h:b] <-(vv_gr_1[i,1:VG$nmi[i]])
    vv_fp[i] <-sum(vv_frac_p1[i,1:VG$nmi[i]])
    vv_fe[i] <-sum(vv_frac_e1[i,1:VG$nmi[i]])
    h=b+1 ; b=b+VG$nmi[i+1]
}

h=1 ;b=VG$nmi[1]
for(i in 1:N_vv){
    vv_dif_2bri[i]<-vv_gr[h+1]-vv_gr[h]
    vv_fp_gr[i]<-vv_fp[i]/VG$nmi[i]
    vv_fe_gr[i]<-vv_fe[i]/VG$nmi[i]
    h=b+1 ; b=b+VG$nmi[i+1]
}

#----------------------------------------------------------------------
# voids R
Nmi_vr=subset(VG$nmi,VG$tipo==0)
igru_vr=subset(VG$igru,VG$tipo==0)

vr_gr_1<-matrix(nrow=N_vr,ncol= 10)
vr_frac_p1<-matrix(nrow=N_vr,ncol= 10)
vr_frac_e1<-matrix(nrow=N_vr,ncol= 10)

for(i in 1:N_vr){
    vr_gr_1[i,1:Nmi_vr[i]]<-subset(HH$mag_r, igru_vr[i]==HH$GId)
    vr_frac_p1[i,1:Nmi_vr[i]]<-subset(HH$passive, igru_vr[i]==HH$GId)
    vr_frac_e1[i,1:Nmi_vr[i]]<-subset(HH$early, igru_vr[i]==HH$GId)
}

vr_gr<-vector('logical',sum(Nmi_vr))
vr_dif_2bri<-vector('logical',N_vr)
vr_fp<-vector('logical',N_vr); vr_fp_gr<-vector('logical',N_vr)
vr_fe<-vector('logical',N_vr); vr_fe_gr<-vector('logical',N_vr)

h=1 ;b=Nmi_vr[1]
for(i in 1:N_vr){
    vr_gr[h:b] <-(vr_gr_1[i,1:Nmi_vr[i]])
    vr_fp[i] <-sum(vr_frac_p1[i,1:Nmi_vr[i]])
    vr_fe[i] <-sum(vr_frac_e1[i,1:Nmi_vr[i]])
    h=b+1 ; b=b+Nmi_vr[i+1]
}

h=1 ;b=Nmi_vr[1]
for(i in 1:N_vr){
    vr_dif_2bri[i]<-vr_gr[h+1]-vr_gr[h]
    vr_fp_gr[i]<-vr_fp[i]/Nmi_vr[i]
    vr_fe_gr[i]<-vr_fe[i]/Nmi_vr[i]
    h=b+1 ; b=b+Nmi_vr[i+1]
}

#----------------------------------------------------------------------
# voids S
Nmi_vs=subset(VG$nmi,VG$tipo==1)
igru_vs=subset(VG$igru,VG$tipo==1)

vs_gr_1<-matrix(nrow=N_vs,ncol= 10)
vs_frac_p1<-matrix(nrow=N_vs,ncol= 10)
vs_frac_e1<-matrix(nrow=N_vs,ncol= 10)

for(i in 1:N_vs){
    vs_gr_1[i,1:Nmi_vs[i]]<-subset(HH$mag_r, igru_vs[i]==HH$GId)
    vs_frac_p1[i,1:Nmi_vs[i]]<-subset(HH$passive, igru_vs[i]==HH$GId)
    vs_frac_e1[i,1:Nmi_vs[i]]<-subset(HH$early, igru_vs[i]==HH$GId)
}

vs_gr<-vector('logical',sum(Nmi_vs))
vs_dif_2bri<-vector('logical',N_vs)
vs_fp<-vector('logical',N_vs); vs_fp_gr<-vector('logical',N_vs)
vs_fe<-vector('logical',N_vs); vs_fe_gr<-vector('logical',N_vs)

h=1 ;b=Nmi_vs[1]
for(i in 1:N_vs){
    vs_gr[h:b] <-(vs_gr_1[i,1:Nmi_vs[i]])
    vs_fp[i] <-sum(vs_frac_p1[i,1:Nmi_vs[i]])
    vs_fe[i] <-sum(vs_frac_e1[i,1:Nmi_vs[i]])
    h=b+1 ; b=b+Nmi_vs[i+1]
}

h=1 ;b=Nmi_vs[1]

for(i in 1:N_vs){
vs_dif_2bri[i]<-vs_gr[h+1]-vs_gr[h]
    vs_fp_gr[i]<-vs_fp[i]/Nmi_vs[i]
    vs_fe_gr[i]<-vs_fe[i]/Nmi_vs[i]
h=b+1 ; b=b+Nmi_vs[i+1]
}

#----------------------------------------------------------------------
# voids IN 
Nmi_vin=subset(VG$nmi,VG$estado==-1)
igru_vin=subset(VG$igru,VG$estado==-1)

vin_gr_1<-matrix(nrow=N_vin,ncol= 10)
vin_frac_p1<-matrix(nrow=N_vin,ncol= 10)
vin_frac_e1<-matrix(nrow=N_vin,ncol= 10)

for(i in 1:N_vin){
    vin_gr_1[i,1:Nmi_vin[i]]<-subset(HH$mag_r, igru_vin[i]==HH$GId)
    vin_frac_p1[i,1:Nmi_vin[i]]<-subset(HH$passive, igru_vin[i]==HH$GId)
    vin_frac_e1[i,1:Nmi_vin[i]]<-subset(HH$early, igru_vin[i]==HH$GId)
}

vin_gr<-vector('logical',sum(Nmi_vin))
vin_dif_2bri<-vector('logical',N_vin)
vin_fp<-vector('logical',N_vin); vin_fp_gr<-vector('logical',N_vin)
vin_fe<-vector('logical',N_vin); vin_fe_gr<-vector('logical',N_vin)

h=1 ;b=Nmi_vin[1]
for(i in 1:N_vin){
    vin_gr[h:b] <-(vin_gr_1[i,1:Nmi_vin[i]])
    vin_fp[i] <-sum(vin_frac_p1[i,1:Nmi_vin[i]])
    vin_fe[i] <-sum(vin_frac_e1[i,1:Nmi_vin[i]])
    h=b+1 ; b=b+Nmi_vin[i+1]
}

h=1 ;b=Nmi_vin[1]
for(i in 1:N_vin){
    vin_dif_2bri[i]<-vin_gr[h+1]-vin_gr[h]
    vin_fp_gr[i]<-vin_fp[i]/Nmi_vin[i]
    vin_fe_gr[i]<-vin_fe[i]/Nmi_vin[i]
    h=b+1 ; b=b+Nmi_vin[i+1]
}

#----------------------------------------------------------------------
# voids EDGE
Nmi_ved=subset(VG$nmi,VG$tipo==1)
igru_ved=subset(VG$igru,VG$tipo==1)

ved_gr_1<-matrix(nrow=N_ved,ncol= 10)
ved_frac_p1<-matrix(nrow=N_ved,ncol= 10)
ved_frac_e1<-matrix(nrow=N_ved,ncol= 10)

for(i in 1:N_ved){
    ved_gr_1[i,1:Nmi_ved[i]]<-subset(HH$mag_r, igru_ved[i]==HH$GId)
    ved_frac_p1[i,1:Nmi_ved[i]]<-subset(HH$passive, igru_ved[i]==HH$GId)
    ved_frac_e1[i,1:Nmi_ved[i]]<-subset(HH$early, igru_ved[i]==HH$GId)
}

ved_gr<-vector('logical',sum(Nmi_ved))
ved_dif_2bri<-vector('logical',N_ved)
ved_fp<-vector('logical',N_ved); ved_fp_gr<-vector('logical',N_ved)
ved_fe<-vector('logical',N_ved); ved_fe_gr<-vector('logical',N_ved)

h=1 ;b=Nmi_ved[1]
for(i in 1:N_ved){
    ved_gr[h:b] <-(ved_gr_1[i,1:Nmi_ved[i]])
    ved_fp[i] <-sum(ved_frac_p1[i,1:Nmi_ved[i]])
    ved_fe[i] <-sum(ved_frac_e1[i,1:Nmi_ved[i]])
    h=b+1 ; b=b+Nmi_ved[i+1]
}

h=1 ;b=Nmi_ved[1]

for(i in 1:N_ved){
ved_dif_2bri[i]<-ved_gr[h+1]-ved_gr[h]
    ved_fp_gr[i]<-ved_fp[i]/Nmi_ved[i]
    ved_fe_gr[i]<-ved_fe[i]/Nmi_ved[i]
h=b+1 ; b=b+Nmi_ved[i+1]
}

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#cairo_pdf("boxplot.pdf")
#pdf("boxplot.pdf")
par(family="serif")
par(cex.lab=1.)       #Tamaño labels
par(cex.main=1.1)        #Tamaño título
par(cex.axis=1.)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(1.8,0.4,0))      	        #margen labels
par(mar=c(4,4,4,4))

par(mfrow = c(6,2))   # dibujo 4x2

par(mar = c(0, 2, 0, 3), oma = c(4, 3, 3, 0)) #pegame este boxplot

#par(mgp=c(1.8,0.4,0))      	        #margen labels




# ------------------------------------------------------
# -----------Sigv----------------------------------
# ------------------------------------------------------
x <- S$sigv #nodos
y <- SS$sigv #grupos FoF
z <- FF$sigv #filamentos
w <- Gf$sigv #campo
v <- VG$sigv #

#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('$\\sigma$ \\[km s$^{-1}$\\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
median(x);median(y);median(z)
# $conf te da los extremos de los notches
boxplot.stats(x, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
boxplot.stats(y, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
# mediana
median(x);median(y);median(z)

# ------------------------------------------------------
# -----------Redshifth----------------------------------
# ------------------------------------------------------
x <-  S$zmedian #nodos
y <- SS$zmedian #grupos FoF
z <- FF$zmedian #filamentos
w <- Gf$zmedian #campo
v <- VG$zmedian #voids

#boxplot(x,y,z,w,v,	
boxplot(x,z,y,w,v,	
	las = 1,notch = TRUE,varwidth = TRUE,outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        #col= c("red",'magenta','darkorange','blue','black'))
        col= c("red",'darkorange','magenta','blue','black'))
mtext(expression(paste(Redshift)), side = 2, cex = 1, line = 3.2, col = "black")
median(x);median(y);median(z)
# $conf te da los extremos de los notches
boxplot.stats(x, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
boxplot.stats(y, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
# mediana
#median(x);median(y);median(z)

# ------------------------------------------------------
# -----------mu----------------------------------
# ------------------------------------------------------
x <-  S$mu  #nodos
y <- SS$mu #grupos FoF
z <- FF$mu #filamentos
w <- Gf$mu #campo
v <- VG$mu #campo
#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= c("red",'darkorange','magenta','blue','black'))

mtext(expression(paste(mu)), side = 2, cex = 1, line = 3.2, col = "black")
median(x);median(y);median(z)
# $conf te da los extremos de los notches
boxplot.stats(x, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
boxplot.stats(y, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
# mediana
median(x);median(y);median(z)

# ------------------------------------------------------
# -----------mabs1----------------------------------
# ------------------------------------------------------
x <- S$rabs1  #nodos
y <- SS$rabs1 #grupos FoF
z <- FF$rabs1 #filamentos
w <- Gf$rabs1 #campo
v <- VG$rabs1 #campo
#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('M_{brigth]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
median(x);median(y);median(z)
# $conf te da los extremos de los notches
boxplot.stats(x, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
boxplot.stats(y, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
# mediana
median(x);median(y);median(z)

# ------------------------------------------------------
# -----------titaG----------------------------------
# ------------------------------------------------------
x <- S$radio_mins #nodos
y <- SS$radio_mins #grupos FoF
z <- FF$radio_mins #filamentos
w <- Gf$radio_mins #campo
v <- VG$radio_mins #voids
#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c('N','GG','F','C'),
        col= c("red",'darkorange','magenta','blue','black'))

laby=TeX('$\\theta$ \\[arcmin\\]') ;
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
median(x);median(y);median(z)
# $conf te da los extremos de los notches
boxplot.stats(x, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
boxplot.stats(y, coef = 1.5, do.conf = TRUE, do.out = FALSE)$conf
# mediana
median(x);median(y);median(z)

# ------------------------------------------------------
# -----------rproy----------------------------------
# ------------------------------------------------------
x <- S$rp  #nodos
y <- SS$rp #grupos FoF
z <- FF$rp #filamentos
w <- Gf$rp #campo
v <- VG$rp #voids
#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c('N','GG','F','C'),
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('$r_p$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")


# ------------------------------------------------------
# -----------dij----------------------------------
# ------------------------------------------------------
x <- S$dij*1000.  #nodos
y <- SS$dij*1000. #grupos FoF
z <- FF$dij*1000. #filamentos
w <- Gf$dij*1000. #campo
v <- VG$dij*1000. #voids

#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------tcr----------------------------------
# ------------------------------------------------------
x <- S$tcr  #nodos
y <- SS$tcr #grupos FoF
z <- FF$tcr #filamentos
w <- Gf$tcr #campo
#w <- (subset(w,w>0.3)
v <- VG$tcr #voids

#boxplot(x,y,z,w,v,
boxplot(x,z,y,w,v,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c('N','GG','F','C','V'),
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('$H_0 \\, t_{cr}$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")




# ------------------------------------------------------
# ----------- M2 - M1 ----------------------------------
# ------------------------------------------------------
boxplot(nodo_dif_2bri,fil_dif_2bri,Fof_dif_2bri,cam_dif_2bri,vv_dif_2bri,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	#names=c(TeX('$Node$'),TeX('$Fil$'),TeX('$Loose$'),TeX('$Non-Emb$'),TeX('$Void$')),
        col= c("red",'darkorange','magenta','blue','black'))

laby=TeX('$m_2 -m_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

x <- S$rabs2 - S$rabs1 #nodos
y <- SS$rabs2 - SS$rabs1 #grupos FoF
z <- FF$rabs2 - FF$rabs1 #filamentos
w <- Gf$rabs2 - Gf$rabs1 #campo
#w <- (subset(w,w>0.3)
v <- VG$rabs2 - VG$rabs1 #voids

boxplot(x,z,y,w,v,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
    col= c("red",'darkorange','magenta','blue','black'))

laby=TeX('$M_2 -M_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")




# ------------------------------------------------------
# ----------- frac passive -----------------------------
# ------------------------------------------------------
boxplot(nodo_fp_gr,fil_fp_gr,Fof_fp_gr,cam_fp_gr,vv_fp_gr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(0,1.1),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$Node$'),TeX('$Fil$'),TeX('$Loose$'),TeX('$Non-Emb$'),TeX('$Void$')),
    col= c("red",'darkorange','magenta','blue','black'))

laby=TeX('$Frac \\, passive \\, glxs$')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# ----------- frac early - -----------------------------
# ------------------------------------------------------

boxplot(nodo_fe_gr,fil_fe_gr,Fof_fe_gr,cam_fe_gr,vv_fe_gr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(0,1.1),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$Node$'),TeX('$Fil$'),TeX('$Loose$'),TeX('$Non-Emb$'),TeX('$Void$')),
        col= c("red",'darkorange','magenta','blue','black'))

laby=TeX('$Frac \\, early \\, glxs$')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")





#dev.off()


#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################

par(ask=TRUE)

#cairo_pdf("boxplot_b.pdf")
par(family="serif")
par(cex.lab=1.)       #Tamaño labels
par(cex.main=1.1)        #Tamaño título
par(cex.axis=1.)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(1.8,0.4,0))      	        #margen labels
par(mar=c(4,4,4,4))

par(mfrow = c(6,2))   # dibujo 4x2

par(mar = c(0, 2, 0, 3), oma = c(4, 3, 3, 0)) #pegame este boxplot

#par(mgp=c(1.8,0.4,0))      	        #margen labels

colores=c("red",'darkorange','magenta','blue','darkgreen','deepskyblue3')


# ------------------------------------------------------
# -----------Sigv----------------------------------
# ------------------------------------------------------
x <- S$sigv #nodos
y <- SS$sigv #grupos FoF
z <- FF$sigv #filamentos
w <- Gf$sigv #campo
v <- VG$sigv #voids
vr <- subset(VG$sigv,VG$tipo==0) #voids
vs <- subset(VG$sigv,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)
laby=TeX('$\\sigma$ \\[km s$^{-1}$\\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------Redshifth----------------------------------
# ------------------------------------------------------
x <-  S$zmedian #nodos
y <- SS$zmedian #grupos FoF
z <- FF$zmedian #filamentos
w <- Gf$zmedian #campo
v <- VG$zmedian #voids
vr <- subset(VG$zmedian,VG$tipo==0) #voids
vs <- subset(VG$zmedian,VG$tipo==1) #voids


boxplot(x,z,y,w,vs,vr,	
	las = 1,notch = TRUE,varwidth = TRUE,outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)
mtext(expression(paste(Redshift)), side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------mu----------------------------------
# ------------------------------------------------------
x <-  S$mu  #nodos
y <- SS$mu #grupos FoF
z <- FF$mu #filamentos
w <- Gf$mu #campo
v <- VG$mu #campo
vr <- subset(VG$mu,VG$tipo==0) #voids
vs <- subset(VG$mu,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)

mtext(expression(paste(mu)), side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------mabs1----------------------------------
# ------------------------------------------------------
x <- S$rabs1  #nodos
y <- SS$rabs1 #grupos FoF
z <- FF$rabs1 #filamentos
w <- Gf$rabs1 #campo
v <- VG$rabs1 #campo
vr <- subset(VG$rabs1,VG$tipo==0) #voids
vs <- subset(VG$rabs1,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)
laby=TeX('M_{brigth]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------titaG----------------------------------
# ------------------------------------------------------
x <- S$radio_mins  #nodos
y <- SS$radio_mins #grupos FoF
z <- FF$radio_mins #filamentos
w <- Gf$radio_mins #campo
v <- VG$radio_mins #voids
vr <- subset(VG$radio_mins,VG$tipo==0) #voids
vs <- subset(VG$radio_mins,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
	#xaxt="n",
        col= colores)

laby=TeX('$\\theta$ \\[arcmin\\]') ;
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------rproy----------------------------------
# ------------------------------------------------------
x <- S$rp  #nodos
y <- SS$rp #grupos FoF
z <- FF$rp #filamentos
w <- Gf$rp #campo
v <- VG$rp #voids
vr <- subset(VG$rp,VG$tipo==0) #voids
vs <- subset(VG$rp,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c('N','GG','F','C'),
        col= colores)
laby=TeX('$r_p$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")


# ------------------------------------------------------
# -----------dij----------------------------------
# ------------------------------------------------------
x <- S$dij*1000.  #nodos
y <- SS$dij*1000. #grupos FoF
z <- FF$dij*1000. #filamentos
w <- Gf$dij*1000. #campo
v <- VG$dij*1000. #voids
vr <- subset(VG$dij*1000.,VG$tipo==0) #voids
vs <- subset(VG$dij*1000.,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
        col= colores)
laby=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------tcr----------------------------------
# ------------------------------------------------------
x <- S$tcr  #nodos
y <- SS$tcr #grupos FoF
z <- FF$tcr #filamentos
w <- Gf$tcr #campo
v <- VG$tcr #voids
vr <- subset(VG$tcr,VG$tipo==0) #voids
vs <- subset(VG$tcr,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)

laby=TeX('$H_0 \\, t_{cr}$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")



# ------------------------------------------------------
# ----------- M2 - M1 ----------------------------------
# ------------------------------------------------------
boxplot(nodo_dif_2bri,fil_dif_2bri,Fof_dif_2bri,cam_dif_2bri,vs_dif_2bri,vr_dif_2bri,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V_S}$'),TeX('$CG_{V_R}$')),
	#names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$VS$'),TeX('$VR$')),
        col= colores)

laby=TeX('$m_2 -m_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

x <-  S$rabs2 -  S$rabs1  #nodos
y <- SS$rabs2 - SS$rabs1  #grupos FoF
z <- FF$rabs2 - FF$rabs1  #filamentos
w <- Gf$rabs2 - Gf$rabs1  #campo
v <- VG$rabs2 - VG$rabs1  #campo
vr <- subset(VG$rabs2-VG$rabs1,VG$tipo==0) #voids
vs <- subset(VG$rabs2-VG$rabs1,VG$tipo==1) #voids

boxplot(x,z,y,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,2.5),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V_S}$'),TeX('$CG_{V_R}$')),
	#names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$VS$'),TeX('$VR$')),
        col= colores)

laby=TeX('$M_2 -M_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")



# ------------------------------------------------------
# ----------- frac passive -----------------------------
# ------------------------------------------------------
boxplot(nodo_fp_gr,fil_fp_gr,Fof_fp_gr,cam_fp_gr,vs_fp_gr,vr_fp_gr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(0,1.1),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$VS$'),TeX('$VR$')),
        col= colores)

laby=TeX('$Frac \\, passive \\, glxs$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# ----------- frac early - -----------------------------
# ------------------------------------------------------

boxplot(nodo_fe_gr,fil_fe_gr,Fof_fe_gr,cam_fe_gr,vs_fp_gr,vr_fp_gr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(0,1.1),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$VS$'),TeX('$VR$')),
        col= colores)

laby=TeX('$Frac \\, early \\, glxs$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")




#dev.off()



#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################

par(ask=TRUE)

#cairo_pdf("boxplot_c.pdf")
par(family="serif")
par(cex.lab=1.)       #Tamaño labels
par(cex.main=1.1)        #Tamaño título
par(cex.axis=1.)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(1.8,0.4,0))      	        #margen labels
par(mar=c(4,4,4,4))

par(mfrow = c(6,2))   # dibujo 4x2

par(mar = c(0, 2, 0, 3), oma = c(4, 3, 3, 0)) #pegame este boxplot

#par(mgp=c(1.8,0.4,0))      	        #margen labels

colores=c("red",'darkorange','magenta','blue','deepskyblue3','darkgreen')


# ------------------------------------------------------
# -----------Sigv----------------------------------
# ------------------------------------------------------
x <- S$sigv #nodos
y <- SS$sigv #grupos FoF
z <- FF$sigv #filamentos
w <- Gf$sigv #campo
v <- VG$sigv #voids
vr <- subset(VG$sigv,VG$estado==0) #voids EDGE
vs <- subset(VG$sigv,VG$estado==-1) #voids IN

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)
laby=TeX('$\\sigma$ \\[km s$^{-1}$\\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------Redshifth----------------------------------
# ------------------------------------------------------
x <-  S$zmedian #nodos
y <- SS$zmedian #grupos FoF
z <- FF$zmedian #filamentos
w <- Gf$zmedian #campo
v <- VG$zmedian #voids
vr <- subset(VG$zmedian,VG$estado==0) #voids
vs <- subset(VG$zmedian,VG$estado==-1) #voids


boxplot(x,z,y,w,vr,vs,	
	las = 1,notch = TRUE,varwidth = TRUE,outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)
mtext(expression(paste(Redshift)), side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------mu----------------------------------
# ------------------------------------------------------
x <-  S$mu  #nodos
y <- SS$mu #grupos FoF
z <- FF$mu #filamentos
w <- Gf$mu #campo
v <- VG$mu #campo
vr <- subset(VG$mu,VG$estado==0) #voids
vs <- subset(VG$mu,VG$estado==-1) #voids

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)

mtext(expression(paste(mu)), side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------mabs1----------------------------------
# ------------------------------------------------------
x <- S$rabs1  #nodos
y <- SS$rabs1 #grupos FoF
z <- FF$rabs1 #filamentos
w <- Gf$rabs1 #campo
v <- VG$rabs1 #campo
vr <- subset(VG$rabs1,VG$estado==0) #voids
vs <- subset(VG$rabs1,VG$estado==-1) #voids

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)
laby=TeX('M_{brigth]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------titaG----------------------------------
# ------------------------------------------------------
x <- S$radio_mins  #nodos
y <- SS$radio_mins #grupos FoF
z <- FF$radio_mins #filamentos
w <- Gf$radio_mins #campo
v <- VG$radio_mins #voids
vr <- subset(VG$radio_mins,VG$estado==0) #voids
vs <- subset(VG$radio_mins,VG$estado==-1) #voids

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
	#xaxt="n",
        col= colores)

laby=TeX('$\\theta$ \\[arcmin\\]') ;
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------rproy----------------------------------
# ------------------------------------------------------
x <- S$rp  #nodos
y <- SS$rp #grupos FoF
z <- FF$rp #filamentos
w <- Gf$rp #campo
v <- VG$rp #voids
vr <- subset(VG$rp,VG$estado==0) #voids
vs <- subset(VG$rp,VG$estado==-1) #voids

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c('N','GG','F','C'),
        col= colores)
laby=TeX('$r_p$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")


# ------------------------------------------------------
# -----------dij----------------------------------
# ------------------------------------------------------
x <- S$dij*1000.  #nodos
y <- SS$dij*1000. #grupos FoF
z <- FF$dij*1000. #filamentos
w <- Gf$dij*1000. #campo
v <- VG$dij*1000. #voids
vr <- subset(VG$dij*1000.,VG$estado==0) #voids
vs <- subset(VG$dij*1000.,VG$estado==-1) #voids

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
        col= colores)
laby=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------tcr----------------------------------
# ------------------------------------------------------
x <- S$tcr  #nodos
y <- SS$tcr #grupos FoF
z <- FF$tcr #filamentos
w <- Gf$tcr #campo
v <- VG$tcr #voids
vr <- subset(VG$tcr,VG$estado==0) #voids
vs <- subset(VG$tcr,VG$estado==-1) #voids

boxplot(x,z,y,w,vr,vs,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
        col= colores)

laby=TeX('$H_0 \\, t_{cr}$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")



# ------------------------------------------------------
# ----------- M2 - M1 ----------------------------------
# ------------------------------------------------------
boxplot(nodo_dif_2bri,fil_dif_2bri,Fof_dif_2bri,cam_dif_2bri,ved_dif_2bri,vin_dif_2bri,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V_{brd}}$'),TeX('$CG_{V_{dnt}}$')),
        col= colores)

laby=TeX('$m_2 -m_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

x <-  S$rabs2 -  S$rabs1  #nodos
y <- SS$rabs2 - SS$rabs1  #grupos FoF
z <- FF$rabs2 - FF$rabs1  #filamentos
w <- Gf$rabs2 - Gf$rabs1  #campo
v <- VG$rabs2 - VG$rabs1  #campo
vr <- subset(VG$rabs2-VG$rabs1,VG$estado==0) #voids
vs <- subset(VG$rabs2-VG$rabs1,VG$estado==-1) #voids

boxplot(x,z,y,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,2.5),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V_{brd}}$'),TeX('$CG_{V_{dnt}}$')),
        col= colores)

laby=TeX('$M_2 -M_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")



# ------------------------------------------------------
# ----------- frac passive -----------------------------
# ------------------------------------------------------
boxplot(nodo_fp_gr,fil_fp_gr,Fof_fp_gr,cam_fp_gr,ved_fp_gr,vin_fp_gr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(0,1.1),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$V_{ED}$'),TeX('$V_{IN}$')),
        col= colores)

laby=TeX('$Frac \\, passive \\, glxs$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# ----------- frac early - -----------------------------
# ------------------------------------------------------

boxplot(nodo_fe_gr,fil_fe_gr,Fof_fe_gr,cam_fe_gr,ved_fe_gr,vin_fe_gr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(0,1.1),
	cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$V_{ED}$'),TeX('$V_{IN}$')),
        col= colores)

laby=TeX('$Frac \\, early \\, glxs$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")



#dev.off()













