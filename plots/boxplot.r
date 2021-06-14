#source("boxplot.r") 
library(bootstrap); library(latex2exp); library(astrolibR)
 #library(magicaxis)
#--------------------------------------------------------------

SS<-read.table("../muestras_finales/compact_in_gg_m3_full")
S<-read.table("../muestras_finales/compact_in_node_m3_full")
Gf<-read.table("../muestras_finales/compact_in_field_m3_full")
FF<-read.table("../muestras_finales/compact_in_filaments_m3_full")
VG<-read.table("../muestras_finales/compact_in_voids_m3_full")

#---------COLNAMES ---------
colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

#R_no<-hist(S$V14, plot=FALSE)
#R_gg<-hist(SS$V11, plot=FALSE) 
#R_fi<-hist(FF$V14, plot=FALSE)
#R_cp<-hist(Gf$V11, plot=FALSE) 

#  plot(R_no$mids,R_no$counts/sum(R_no$counts),type='s',col='red',lwd=4)
#points(R_gg$mids,R_gg$counts/sum(R_gg$counts),type='s',col='magenta',lwd=4)
#points(R_fi$mids,R_fi$counts/sum(R_fi$counts),type='s',col='orange',lwd=4)
#points(R_cp$mids,R_cp$counts/sum(R_cp$counts),type='s',col='blue',lwd=4)

#pdf("boxplot.pdf")
par(family="serif")
par(cex.lab=1.)       #Tamaño labels
par(cex.main=1.1)        #Tamaño título
par(cex.axis=1.)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(1.8,0.4,0))      	        #margen labels
par(mar=c(4,4,4,4))

par(mfrow = c(4,2))   # dibujo 4x2

par(mar = c(0, 2, 0, 3), oma = c(4, 3, 3, 0)) #pegame este boxplot

# ------------------------------------------------------
# -----------Redshifth----------------------------------
# ------------------------------------------------------
z_max=0.087
x <-  subset(S$zmedian,S$zmedian <= z_max) #nodos
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
# -----------Sigv----------------------------------
# ------------------------------------------------------
x <- subset(S$sigv, S$zmedian <= z_max) #nodos
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
# -----------mu----------------------------------
# ------------------------------------------------------
x <-  subset(S$mu, S$zmedian <= z_max)  #nodos
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
x <- subset(S$rabs1, S$zmedian <= z_max)  #nodos
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
x <- subset(S$radio_mins, S$zmedian <= z_max)  #nodos
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
x <- subset(S$rp, S$zmedian <= z_max)  #nodos
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
x <- subset(S$dij*1000., S$zmedian <= z_max)  #nodos
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
	names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('$d_{ij}$ \\[kpc h$^{-1}$ \\]')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")

# ------------------------------------------------------
# -----------tcr----------------------------------
# ------------------------------------------------------
x <- subset(S$tcr, S$zmedian <= z_max)  #nodos
y <- SS$tcr #grupos FoF
z <- FF$tcr #filamentos
w <- Gf$tcr #campo
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
	#xaxt="n",
	#names=c('N','GG','F','C','V'),
	names=c(TeX('$CG_{N}$'),TeX('$CG_{F}$'),TeX('$CG_{L}$'),TeX('$CG_{C}$'),TeX('$CG_{V}$')),
        col= c("red",'darkorange','magenta','blue','black'))
laby=TeX('$H_0 \\, t_{cr}$ ')
mtext(laby, side = 2, cex = 1, line = 3.2, col = "black")
#dev.off()














