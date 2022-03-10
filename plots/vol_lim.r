#source("vol_lim.r")

library(latex2exp)
library(magicaxis)

SS<-read.table("../data/compact_in_gg_m3_full")
S<-read.table("../data/compact_in_node_m3_full")
Gf<-read.table("../data/compact_in_field_m3_full")
FF<-read.table("../data/compact_in_filaments_m3_full")
VG<-read.table("../data/compact_in_voids_m3_full")

#HH<-read.table("../catalogos/tab_gal_gru.dat")
#gal<-read.table("catalogos/DR12_tempel.dat")
gal<-read.table("../../galgrufof_z_Mr.dat") # vol_samples_new.f


colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

#colnames(HH)[c(1,2,3,4,5,6,7)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g")

z_no<-S$zmedian; R_no<-S$rabs1 

z_gg<-SS$zmedian; R_gg<-SS$rabs1

z_fi<-FF$zmedian; R_fi<-FF$rabs1

z_cp<-Gf$zmedian; R_cp<-Gf$rabs1

z_vv<-VG$zmedian; R_vv<-VG$rabs1

library(dplyr)
df=data.frame( gal$V1, gal$V3, gal$V4, gal$V5, gal$V6)
new_df=sample_frac(df,0.2)

dim(new_df); nrow(new_df); ncol(new_df)

zgal <- new_df[,1]
Mgal <- new_df[,2]
mag_env <- new_df[,3]
mag_env_cg=mag_env-3

Mgal_kcorr <- new_df[,4]
mag_env_kcorr <- new_df[,5]
mag_env_kcorr_cg=mag_env_kcorr-3

#pdf("vol_lim.pdf")


par(mar=c(5,5,1,1))  #c(b,l,t,r)
#par(oma=c(0,0,0,0))  #c(b,l,t,r)

plot(zgal,Mgal,ylim=c(-17,-23.2),xlim=c(0,0.19),pch='.',cex=0.5,col='gray',
     xlab='Redshift',ylab=TeX('$M_{bri}-5log_{10}(h)$'),cex.lab=1.8, cex.axis=1.5)
points(z_no,R_no,pch=16,cex=0.6,col='red')
points(z_fi,R_fi,pch=5,cex=0.6,col='darkorange')
points(z_gg,R_gg,pch=0,cex=0.6,col='magenta')
points(z_cp,R_cp,pch=4,cex=0.6,col='blue')
points(z_vv,R_vv,pch=11,cex=0.6,col='black')
magaxis(label=FALSE)

abline(v=0.01,lty=2)
#abline(v=0.1)
#abline(h=-19.77)

segments(-0.5, -19.77, 0.1 , -19.77,col='black',lty=2)
segments(0.1, -19.77, 0.1 , -23.77,col='black',lty=2)

#envolvente
zt=c(0.01,0.013,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.195,0.2)
med_zz=numeric()
med_mm=numeric()
for (i in 1:length(zt)-1){
	print(i)
zz=subset(zgal,zgal>zt[i] & zgal<zt[i+1])
mm=subset(mag_env,zgal>zt[i] & zgal<zt[i+1])
med_zz[i]=median(zz)
med_mm[i]=median(mm)
}

points(med_zz,med_mm, type="l",pch='.')
lines(med_zz,med_mm, type="l", col="black", lty=1)
lines(med_zz,med_mm-3, type="l", col="black", lty=1)


legend(0.12,-19.7, legend=c("Nodes","Filament", "Loose Groups","Non-Embedded","Voids"),col=c('red',"darkorange",'magenta',"blue","black"), pch=c(16,5,0,4,11), cex=1.4,border = NULL,bty='n')

#dev.off()


#z_v1=0.1
#rrglx=292.678683
#dl=rrglx*(1.+z_v1)
#M_v1=17.77-5.*log10(dl)-25

