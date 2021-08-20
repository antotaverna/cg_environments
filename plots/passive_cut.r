library(bootstrap); library(latex2exp); library(astrolibR)
#library(magicaxis)
library(cosmoFns)

G_FF<-read.table("../../catalogos/galgrufof_408.dat")
gal_fof<-read.table("../../catalogos/mabs_galgrufof.dat")
tempel<-read.table("../../tempel_z_Mr.dat") #empieza en z=0.1 guarda!!

colnames(G_FF)[c(1,2,3,4,5,6,7)] <- c("igru",'ra','dec','z','r','g','iggrup')
z_apr<-G_FF$z
r_apr<-G_FF$r
g_apr<-G_FF$g


#M_r1<-gal_fof$V8
#M_g1<-gal_fof$V9
M_r1<-tempel$V4
M_g1<-tempel$V5
M_r[1:length(z_apr)] =0
M_g[1:length(z_apr)] =0


M_r<-subset(M_r1,M_r1 >=-23. & M_r1 <= -17.)
M_g<-subset(M_g1,M_r1 >=-23. & M_r1 <= -17.)

color1= M_g-M_r
col_a= subset(color1,M_r1 >=-23. & M_r1 <= -22.)
col_b= subset(color1,M_r1 >=-22. & M_r1 <= -21.)
col_c= subset(color1,M_r1 >=-21. & M_r1 <= -20.)
col_c= subset(color1,M_r1 >=-20. & M_r1 <= -18.)


#------------------------------------------------------------------------------------------------------
G<-read.table("../data/tab_gal_gru.dat")
colnames(G)[c(1,2,3,4,5,6,7,8,9)
            ] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g","Mag_r","Mag_g")
color=G$Mag_g- G$Mag_r

#------------------------------------------------------------------------------------------------------
#   PLOT    --------------
#------------------------------------------------------------------------------------------------------
#png("col_tempel.png")
library(dplyr)
par(mfrow=c(2,3))
par(mar=c(5,3,1,1))  #c(b,l,t,r)
par(oma=c(0,1,1,1))  #c(b,l,t,r)

limx=c(0.2,1.1)


hist(color1,breaks=seq(from = -3.0, to = 7.1, by = 0.025),xlim=limx,main='')
hist(col_a,breaks=seq(from = -3.0, to = 7.1, by = 0.025),xlim=limx,main='')
legend(0.3,200,legend=c('[-23,-22]'),col='green')
hist(col_b,breaks=seq(from = -3.0, to = 7.1, by = 0.025),xlim=limx,main='')
legend(0.3,5000,legend=c('[-22,-21]'),col='green')
hist(col_c,breaks=seq(from = -3.0, to = 7.1, by = 0.025),xlim=limx,main='')
legend(0.3,6000,legend=c('[-21,-20]'),col='green')
hist(col_d,breaks=seq(from = -3.0, to = 7.1, by = 0.05),xlim=limx,main='')
legend(0.3,10000,legend=c('[-20,-18]'),col='green')
#hist(col_e,breaks=seq(from = -3.0, to = 7.1, by = 0.05),xlim=limx,main='')
#legend(0.3,400,legend=c('[-19,-18]'),col='green')
#hist(col_f,breaks=seq(from = -3.0, to = 7.1, by = 0.05),xlim=limx,main='')
#legend(0.3,200,legend=c('[-18,-17]'),col='green')





M_tol<-rbind(c(M_g,G$Mag_r))
col_tol<-rbind(c(color1 ,color ))

df=data.frame(M_g,color1)
new_df=sample_frac(df,0.3)
M_g_n <- new_df[,1]
color1_n <- new_df[,2]


labx=TeX(' $\\M_r-5 \\, \\log{(h)}}\\]$')
plot(NA,xlim=c(-23,-17.5),ylim=c(0.2,1.1),xlab=labx,ylab='g-r',main='')
points(M_g_n,color1_n,pch=16,cex=0.1,col='grey51')
#points(G$Mag_r,color,pch=16,cex=0.2,col='red')

mycol <- rgb(1, 0, 0, alpha=0.5) # 50% transparent black
require(KernSmooth)

l<-data.frame(G$Mag_r,color)
colnames(l)[c(1,2)] <-c("Mag", "col")

est <- bkde2D(l[c("Mag", "col")], bandwidth=c(0.05, 0.05), gridsize=c(100, 100))
#--Display as a contour map:
#with(est, contour(x1, x2, fhat, drawlabels=TRUE, add=TRUE, col=mycol,labcex=1))
with(est, contour(x1, x2, fhat, drawlabels=TRUE, add=TRUE, col="black",labcex=1))

#y=0.73-0.03*(G$Mag_r+20)
#points(G$Mag_r,y,type='l',col='black')
#y=0.72-0.03*(G$Mag_r+20)
#points(G$Mag_r,y,type='l',col='green')
y=0.7-0.03*(G$Mag_r+20)
points(G$Mag_r,y,type='l',col='blue')
y=0.8-0.03*(G$Mag_r+20)
points(G$Mag_r,y,type='l',col='magenta')
abline(h=0.7,lw=1, lty=1,col='darkgreen')

legend(-23.0, 0.3, legend=c("Blanton+07","g-r = 0.7","g-r=0.7-0.03*(Mr+20)"),
       col=c("magenta","darkgreen", "blue"), lty=1, cex=0.8)




#dev.off()



#png("col_tempel.png")
M_r1<-tempel$V4
M_g1<-tempel$V5

par(mfrow=c(2,1))
limx=c(0.2,1.5)
hist(M_g1-M_r1,breaks=1000,xlim=limx,main='')
hist(tempel$V3-tempel$V2,breaks=1000,xlim=limx,main='')

