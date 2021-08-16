library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
library(cosmoFns)

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


G_FF<-read.table("galgrufof_408.dat")

colnames(G_FF)[c(1,2,3,4,5,6,7)] <- c("igru",'ra','dec','z','r','g','iggrup')
z_apr<-G_FF$z
r_apr<-G_FF$r
g_apr<-G_FF$g
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
M_r1<-vector('logical',length(z_apr)) ; M_r[1:length(z_apr)] =0
M_g1<-vector('logical',length(z_apr)) ; M_g[1:length(z_apr)] =0


for(l in 1:length(z_apr)){
 dL<-D.L(z_apr[l], omega.m =0.315, omega.lambda =0.67, H.0 = 67.3)

 M_r1[l]= r_apr[l]-25.-5.*log10(dL)-K_r2[l]
 M_r1[l]=M_r1[l]-5*log10(h_l)
 M_g1[l]= g_apr[l]-25.-5.*log10(dL)-K_g2[l]
 M_g1[l]=M_g1[l]-5*log10(h_l)
}

M_g<-subset(M_g1,M_r1 >=-23. & M_r1 <= -17. & M_g1 >=-23. & M_g1 <= -17.)
M_r<-subset(M_r1,M_r1 >=-23. & M_r1 <= -17. & M_g1 >=-23. & M_g1 <= -17.)

color1= M_g-M_r

G<-read.table("cg_environments/data/tab_gal_gru.dat")
colnames(G)[c(1,2,3,4,5,6,7,8,9)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g","Mag_r","Mag_g")
color=G$Mag_g- G$Mag_r

#------------------------------------------------------------------------------------------------------
library(dplyr)

#pdf("passive_cut.pdf")
M_tol<-rbind(c(M_g,G$Mag_r))
col_tol<-rbind(c(color1 ,color ))

df=data.frame(M_g,color1)
new_df=sample_frac(df,0.1)
M_g_n <- new_df[,1]
color1_n <- new_df[,2]

labx=TeX(' $\\M_r-5 \\, \\log{(h)}}\\]$')
#plot(NA,xaxt = "n",axes=FALSE,xlim=c(-17,-23),ylim=c(0.2,1.1),xlab=labx,ylab='g-r')
plot(NA,xlim=c(-23,-17.5),ylim=c(0.2,1.1),xlab=labx,ylab='g-r')
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

