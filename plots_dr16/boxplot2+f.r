library(bootstrap); library(latex2exp); library(astrolibR)
 library(magicaxis)
#--------------------------------------------------------------

HH<-read.table("../../structures/enviros_dr16/glx_CG_members+faints.dat") #flux limited


SS<-read.table("../../structures/enviros_dr16/compact_in_loose_v1+f.dat")
S<-read.table("../../structures/enviros_dr16/compact_in_nodes_v1+f.dat")
Gf<-read.table("../../structures/enviros_dr16/compact_in_NE_v1+f.dat")
FF<-read.table("../../structures/enviros_dr16/compact_in_filaments_v1+f.dat")
VR<-read.table("../../structures/enviros_dr16/compact_in_voidsR_v1+f.dat")
VS<-read.table("../../structures/enviros_dr16/compact_in_voidsS_v1+f.dat")


#---------COLNAMES ---------

colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)] <- c('igru','alcm','delcm','zmedian','nmi','radio_mins','mu','rmag_bri','sigv','rp','mvir','rabs1','rabs2','lgroup','tcr','dij','flag','fr_red','fr_early')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)] <- c("igru",'alcm','delcm','zmedian','nmi','radio_mins','mu','rmag_bri','sigv','rp','mvir','rabs1','rabs2','lgroup','tcr','dij','flag','fr_red','fr_early')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)] <- c("igru",'alcm','delcm','zmedian','nmi','radio_mins','mu','rmag_bri','sigv','rp','mvir','rabs1','rabs2','lgroup','tcr','dij','flag','fr_red','fr_early')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)] <- c("igru",'alcm','delcm','zmedian','nmi','radio_mins','mu','rmag_bri','sigv','rp','mvir','rabs1','rabs2','lgroup','tcr','dij','flag','fr_red','fr_early')

colnames(VR)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)] <- c("igru",'alcm','delcm','zmedian','nmi','radio_mins','mu','rmag_bri','sigv','rp','mvir','rabs1','rabs2','lgroup','tcr','dij','flag','fr_red','fr_early')
colnames(VS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)] <- c("igru",'alcm','delcm','zmedian','nmi','radio_mins','mu','rmag_bri','sigv','rp','mvir','rabs1','rabs2','lgroup','tcr','dij','flag','fr_red','fr_early')


#colnames(HH)[c(1,2,3,4,5,6,7,8,9)] <-c("galid", "RA", "Dec", "mag_r", "Redshift", "rabs","color","cindex", "GId", "passive", "early")


#------------------------------------------------------------------
#------------------------------------------------------------------
#seleccion de galaxias miembros
#------------------------------------------------------------------

N_nodo=length(S$igru)
N_Fof=length(SS$igru)
N_fi=length(FF$igru)
N_cam=length(Gf$igru)
N_vr=length(VR$igru)
N_vs=length(VS$igru)
N_vv=N_vr+N_vs


#cairo_pdf("figuras/boxplot2+f.pdf")
#par(family="serif")
par(cex.lab=1.9)       #Tamaño labels
par(cex.axis=1.2)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(4,0.4,0))  #margen labels
par(mar=c(0,0,0,0))
par(mfrow = c(4,2))   # dibujo 4x2

par(mar = c(0, 2.5, 0, 3.5), oma = c(4, 3, 3, 0)) #pegame este boxplot

#par(mgp=c(1.8,0.4,0))      	        #margen labels

colores=c("red",'darkorange','magenta','blue','darkgreen','deepskyblue3')


#--------- prueba color
a=c("red2",'orange2','magenta3','blue3','green4','deepskyblue3')
a=c("#F8766D","#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
a1 <- col2rgb(a)

#alpha=0.5
#cn =rgb(a1[1,1],a1[2,1],a1[3,1], alpha=alpha, maxColorValue = 255)
#clg=rgb(a1[1,2],a1[2,2],a1[3,2], alpha=alpha, maxColorValue = 255)
#cf =rgb(a1[1,3],a1[2,3],a1[3,3], alpha=alpha, maxColorValue = 255)
#cne=rgb(a1[1,4],a1[2,4],a1[3,4], alpha=alpha, maxColorValue = 255)
#cvs=rgb(a1[1,5],a1[2,5],a1[3,5], alpha=alpha, maxColorValue = 255)
#cvr=rgb(a1[1,6],a1[2,6],a1[3,6], alpha=alpha, maxColorValue = 255)

#col_pru=c(cn,clg,cf,cne,cvs,cvr)

# transform to HSV space
a2 <- rgb2hsv(a1)
# you can try different scaling values e.g. between 0.3 - 0.6
n <- 0.45
col_pru=hsv(a2[1,], a2[2,]*n, a2[3,])
#---------------------


# ------------------------------------------------------
# -----------Sigv----------------------------------
# ------------------------------------------------------
x <- S$sigv #nodos
y <- SS$sigv #grupos FoF
z <- FF$sigv #filamentos
w <- Gf$sigv #campo
#v <- VG$sigv #voids
vr <- VR$sigv #voids
vs <- VS$sigv #voids

#boxplot(log10(x),log10(z),log10(y),log10(w),log10(vs),log10(vr),
boxplot(x, y, z, w, vs, vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	whisklty = 0, staplelty=0,
	outline = FALSE,
    whiskers = FALSE,
    ylim = c(100,600),
    #cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n", yaxt="n",
    log='y',
    col= col_pru)
    #col= colores)

laby=TeX('$\\sigma_v$ \\[km s$^{-1}$\\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
rango =c(2, 2.3, 2.6, 2.8)
#mtext(c(TeX('$100$'),TeX('$200$'),TeX('$400 $'),TeX('$ 700 $')),at=rango, side =2, cex = 1, line = 1.2, col = "black")



magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)#,at = rango)
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)

# ------------------------------------------------------
# -----------dij----------------------------------
# ------------------------------------------------------
x <- S$dij  #nodos
y <- SS$dij #grupos FoF
z <- FF$dij #filamentos
w <- Gf$dij #campo
#v <- VG$dij #voids
vr <- VR$dij #voids
vs <- VS$dij #voids

#boxplot(log10(x),log10(z),log10(y),log10(w),log10(vs),log10(vr),	
boxplot(x,y,z,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,
    ylim = c(35,160),
    whisklty = 0, staplelty=0,
#	lwd=c(1.2,1.2,1.5),
    #cex.axis=1.2,
    log = "y",
	xaxt="n", yaxt="n",
    col= col_pru)

laby=TeX('$<d_{ij}>$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

rango =c( 100)

#mtext(c(TeX('$50$'),TeX('$80$'),TeX('$100 $'),TeX('$ 150 $')),at=rango, side =2, cex = 1, line = 1.2, col = "black")

#magaxis(side=2, las=2,at = rango)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2, logpretty=TRUE) #, at = rango)
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)


# ------------------------------------------------------
# -----------mabs1----------------------------------
# ------------------------------------------------------
x <- S$rabs1  #nodos
y <- SS$rabs1 #grupos FoF
z <- FF$rabs1 #filamentos
w <- Gf$rabs1 #campo
#v <- VG$rabs1 #campo
vr <- VR$rabs1
vs <- VS$rabs1 #voids

boxplot(x,y,z,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,	whisklty = 0, staplelty=0,
    ylim = c(-22.6,-20.6),
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",yaxt="n",
    col= col_pru)
laby=TeX('$M_{bri} $')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
magaxis(side=2, majorn=4, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)


# ------------------------------------------------------
# -----------mu----------------------------------
# ------------------------------------------------------
x <-  S$mu  #nodos
y <- SS$mu #grupos FoF
z <- FF$mu #filamentos
w <- Gf$mu #campo
#v <- VG$mu #campo
vr <- VR$mu
vs <- VS$mu #voids

boxplot(x,y,z,w,vs,vr,	
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,	whisklty = 0, staplelty=0,
    ylim = c(22.5,26.3),
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",yaxt="n",
    col= col_pru)

mtext(TeX('$\\mu$\\[mag  $arcsec^{-2}\\]$'), side = 2, line = 3.2, col = "black")
magaxis(side=2, majorn=3, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)

# ------------------------------------------------------
# ----------- M2 - M1 ----------------------------------
# ------------------------------------------------------

x <-  S$rabs2 -  S$rabs1  #nodos
y <- SS$rabs2 - SS$rabs1  #grupos FoF
z <- FF$rabs2 - FF$rabs1  #filamentos
w <- Gf$rabs2 - Gf$rabs1  #campo
#v <- VG$rabs2 - VG$rabs1  #campo
vr <- VR$rabs2 -  VR$rabs1
vs <- VS$rabs2 -  VS$rabs1

boxplot(x,y,z,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,	whisklty = 0, staplelty=0,
    ylim = c(0.3,2.2),
	#cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",yaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V_S}$'),TeX('$CG_{V_R}$')),
	#names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$VS$'),TeX('$VR$')),
        col= col_pru)

laby=TeX('$M_2$ - $M_1$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)

# ------------------------------------------------------
# -----------tcr----------------------------------
# ------------------------------------------------------
x <- S$tcr  #nodos
y <- SS$tcr #grupos FoF
z <- FF$tcr #filamentos
w <- Gf$tcr #campo
#v <- VG$tcr #voids
vr <- VR$tcr
vs <- VS$tcr

#boxplot(log10(x),log10(y),log10(z),log10(w),log10(vs),log10(vr),	
boxplot(x,y,z,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,
	outline = FALSE,	whisklty = 0, staplelty=0,
    ylim = c(0.007,0.1),
#	lwd=c(1.2,1.2,1.5),
	xaxt="n", yaxt="n",
    log='y',
    col= col_pru)

laby=TeX('$H_0 \\, t_{cr}$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

#rango =c(-2, -1.69, -1.39, -1.15)
#mtext(c(TeX('$0.01$'),TeX('$0.02$'),TeX('$0.04 $'),TeX('$ 0.07 $')),at=rango, side =2, cex = 1, line = 1.2, col = "black")

magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2,at = rango)
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)

# ------------------------------------------------------
# ----------- frac passive -----------------------------
# ------------------------------------------------------
x <- S$fr_red  #nodos
y <- SS$fr_red #grupos FoF
z <- FF$fr_red #filamentos
w <- Gf$fr_red #campo
#v <- VG$fr_red #voids
vr <- VR$fr_red
vs <- VS$fr_red

boxplot(x,y,z,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,	whisklty = 0, staplelty=0,
	outline = FALSE,
    ylim = c(0.2,1.1),
	#cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",yaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$SV$'),TeX('$RV$')),
    col= col_pru)

laby=TeX('$F_{red}$ ')
mtext(laby, side = 2, cex = 1.3, line = 3.2, col = "black")
mtext(c(TeX('$CG_N$'),TeX('$CG_{LG}$'),TeX('$CG_{F}$'),TeX('$ CG_{NE} $'),TeX('$CG_{VS} $'),TeX('$CG_{VR}$')),
    at=c(1,2,3,4,5,6), side =1, cex = 0.9, line = 1.2, col = "black")
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)


# ------------------------------------------------------
# ----------- frac early - -----------------------------
# ------------------------------------------------------
x <- S$fr_early  #nodos
y <- SS$fr_early #grupos FoF
z <- FF$fr_early #filamentos
w <- Gf$fr_early #campo
#v <- VG$fr_early #voids
vr <- VR$fr_early
vs <- VS$fr_early

boxplot(x,y,z,w,vs,vr,
	las = 1,
	notch = TRUE,
	varwidth = TRUE,	whisklty = 0, staplelty=0,
	outline = FALSE,
    ylim = c(0.2,1.1),
	#cex.axis=1.2,
#	lwd=c(1.2,1.2,1.5),
	xaxt="n",yaxt="n",
	#names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
	#names=c(TeX('$N$'),TeX('$F$'),TeX('$L$'),TeX('$NE$'),TeX('$VS$'),TeX('$VR$')),
        col= col_pru)

laby=TeX('$F_{early}$ ')
mtext(laby, side = 2, cex = 1.3, line = 3.2, col = "black")
mtext(c(TeX('$CG_N$'),TeX('$CG_{LG}$'),TeX('$CG_{F}$'),TeX('$ CG_{NE} $'), TeX('$CG_{VS} $'),TeX('$CG_{VR}$')),
    at=c(1,2,3,4,5,6), side =1, cex = 0.9, line = 1.2, col = "black")
magaxis(side=1, majorn=5, minorn=5, tcl=0.6, ratio=0., labels=FALSE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)

#dev.off()
