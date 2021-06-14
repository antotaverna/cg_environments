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
fil_gr_1[i,1:FF$nmi[i]]<-subset(HH$mag_r, FF$igru[i]==HH$GId)
}

fil_gr<-vector('logical',sum(FF$nmi)); fil_dif_2bri<-vector('logical',N_fi)

h=1 ;b=FF$nmi[1]

for(i in 1:N_fi){
fil_gr[h:b] <-(fil_gr_1[i,1:FF$nmi[i]])
h=b+1 ; b=b+FF$nmi[i+1]

}
h=1 ;b=FF$nmi[1]

for(i in 1:N_fi){
fil_dif_2bri[i]<-fil_gr[h+1]-fil_gr[h]
h=b+1 ; b=b+FF$nmi[i+1]
}

#----------------------------------------------------------------------

# Nodos
nodo_gr_1<-matrix(nrow=N_nodo,ncol= 10)
for(i in 1:N_nodo){
nodo_gr_1[i,1:S$nmi[i]]<-subset(HH$mag_r, S$igru[i]==HH$GId)
}

nodo_gr<-vector('logical',sum(S$nmi)); nodo_dif_2bri<-vector('logical',N_nodo)
h=1 ;b=S$nmi[1]

for(i in 1:N_nodo){
nodo_gr[h:b] <-(nodo_gr_1[i,1:S$nmi[i]])
h=b+1 ; b=b+S$nmi[i+1]
}

h=1 ;b=S$nmi[1]
for(i in 1:N_nodo){
nodo_dif_2bri[i]<-nodo_gr[h+1]-nodo_gr[h]
h=b+1 ; b=b+S$nmi[i+1]
}

#----------------------------------------------------------------------
# FOF
Fof_gr_1<-matrix(nrow=N_Fof,ncol= 10)
for(i in 1:N_Fof){
Fof_gr_1[i,1:SS$nmi[i]]<-subset(HH$mag_r, SS$igru[i]==HH$GId)
}

Fof_gr<-vector('logical',sum(SS$nmi)); Fof_dif_2bri<-vector('logical',N_Fof)
h=1 ;b=SS$nmi[1]

for(i in 1:N_Fof){
Fof_gr[h:b] <-(Fof_gr_1[i,1:SS$nmi[i]])
h=b+1 ; b=b+SS$nmi[i+1]
}
h=1 ;b=SS$nmi[1]

for(i in 1:N_Fof){
Fof_dif_2bri[i]<-Fof_gr[h+1]-Fof_gr[h]
h=b+1 ; b=b+SS$nmi[i+1]
}


#----------------------------------------------------------------------
# Campo
cam_gr_1<-matrix(nrow=N_cam,ncol= 10)
for(i in 1:N_cam){
cam_gr_1[i,1:Gf$nmi[i]]<-subset(HH$mag_r, Gf$igru[i]==HH$GId)
}

cam_gr<-vector('logical',sum(Gf$nmi)); cam_dif_2bri<-vector('logical',N_cam)
h=1 ;b=Gf$nmi[1]

for(i in 1:N_cam){
cam_gr[h:b] <-(cam_gr_1[i,1:Gf$nmi[i]])
h=b+1 ; b=b+Gf$nmi[i+1]
}

h=1 ;b=Gf$nmi[1]

for(i in 1:N_cam){
cam_dif_2bri[i]<-cam_gr[h+1]-cam_gr[h]
h=b+1 ; b=b+Gf$nmi[i+1]
}


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
# voids
vv_gr_1<-matrix(nrow=N_vv,ncol= 10)
for(i in 1:N_vv){
vv_gr_1[i,1:VG$nmi[i]]<-subset(HH$mag_r, VG$igru[i]==HH$GId)
}

vv_gr<-vector('logical',sum(VG$nmi)); vv_dif_2bri<-vector('logical',N_vv)
h=1 ;b=VG$nmi[1]

for(i in 1:N_vv){
vv_gr[h:b] <-(vv_gr_1[i,1:VG$nmi[i]])
h=b+1 ; b=b+VG$nmi[i+1]
}

h=1 ;b=VG$nmi[1]

for(i in 1:N_vv){
vv_dif_2bri[i]<-vv_gr[h+1]-vv_gr[h]
h=b+1 ; b=b+VG$nmi[i+1]
}


#------------------------------------------------------------------------------------------------------------


v1 <- c('Voids','Field','Fil.','FOF','Node')
v2 <- c(1:5) ; v3 <- c(1.1:4.1) ;v4 <- c(1.2:4.2)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
laby=TeX('$\\sigma$') 
labx=TeX('$M_r$')




par(family="serif")
par(cex.lab=1.2)       #Tamaño labels
par(cex.main=1.2)        #Tamaño título
par(cex.axis=1.)      #Tamaño números en ejes.
par(lwd=2)
par(cex=2)
par(mgp=c(1.8,0.4,0))      	        #margen labels
par(mar=c(4,4,4,4))
#-------------------------------------------------------------
rojo    <- rgb(1,0,0,1)
naranja <- rgb(1, 101/255, 0,1)
verde   <- rgb(24/255, 98/255, 24/255,1)
azul    <- rgb(0,0,1,1)
#--------------------------------------------------------------

boxplot(fil_dif_2bri,nodo_dif_2bri,Fof_dif_2bri,cam_dif_2bri,vv_dif_2bri,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
        #ylim = c(0,0.2),
	cex.axis=0.8,
#	lwd=c(1.2,1.2,1.5),
	#xaxt="n",
	names=c('Node','FoF','Fil','Field','Voids'),
        col= c("red",'magenta','darkorange','blue','black'))
laby=TeX('$M_2 -M_1$ ')
mtext(laby, side = 2, cex = 1.5, line = 3.2, col = "black")








