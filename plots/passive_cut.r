library(bootstrap); library(latex2exp); library(astrolibR)
#library(magicaxis)
library(cosmoFns)

tempel<-read.table("../data/tempel_Mr.dat")

M_r1<-tempel$V7
M_g1<-tempel$V8

M_r<-subset(M_r1,M_r1 >=-23. & M_r1 <= -17.)
M_g<-subset(M_g1,M_r1 >=-23. & M_r1 <= -17.)

color1= M_g-M_r

#------------------------------------------------------------------------------------------------------
#   PLOT    --------------
#------------------------------------------------------------------------------------------------------
#pdf("passive_cut.pdf")
library(dplyr)
library(TeX)
par(mfrow=c(1,1))
par(mar=c(5,3,1,1))  #c(b,l,t,r)
par(oma=c(0,1,1,1))  #c(b,l,t,r)

#hist(color1,breaks=seq(from = -3.0, to = 7.1, by = 0.025),xlim=c(0.2,1.1),main='')

df=data.frame(M_r,color1)
new_df=sample_frac(df,0.05)
M_g_n <- new_df[,1]
color1_n <- new_df[,2]

#labx=TeX(' $\\M_r-5 \\, \\log{(h)}}\\]$')
labx="M_r-5log(h)"
plot(NA,xlim=c(-23,-18.5),ylim=c(0.3,1.0),xlab=labx,ylab='g-r',main='')
points(M_g_n,color1_n,pch=16,cex=0.1,col='grey40')
#points(G$Mag_r,color,pch=16,cex=0.2,col='red')

#---------------------
##--Display as a contour map:
#---------------------
#require(KernSmooth)
#l<-data.frame(M_r,color1)
#colnames(l)[c(1,2)] <-c("Mag", "col")
#est <- bkde2D(l[c("Mag", "col")], bandwidth=c(0.05, 0.05), gridsize=c(100, 100))
#with(est, contour(x1, x2, fhat, drawlabels=TRUE, add=TRUE, col="black",labcex=1))


#Ajustes Lineales
#y=0.7-0.03*(G$Mag_r+20)
#points(G$Mag_r,y,type='l',col='blue')
y=0.8-0.03*(G$Mag_r+20)
points(G$Mag_r,y,type='l',col='magenta')
abline(h=0.740871152583451,lw=1, lty=1,col='darkgreen')

nn=read.table("gmm_intersect.dat")
xx=c(-22.5,-21.75,-21.25,-20.75,-20.25,-19.75,-19.25)#,-18.5)
yy=c(nn$V2[2],nn$V2[3],nn$V2[4],nn$V2[5],nn$V2[6],nn$V2[7],nn$V2[8])#,nn$V2[9])
points(xx,yy,col="blue",pch=16)
#points(xx,yy,col="blue",type="l")


#Ajuste cuadratico
p=poly(xx,yy, degree = 2)
modelo_p <- lm(p)

model <- lm(yy ~ poly(xx,2))
coef1 <- lm(yy ~ xx + I(xx^2))$coef
lines(xseq,  coef1[3]*xseq^2+ coef1[2]*xseq+ coef1[1], col='red')


legend(-23.0, 0.45, legend=c("Blanton+07","Total range [-23,-17]","Nuestro spliteador","fit cuadratico"),
       col=c("magenta","darkgreen", "blue","red"), lty=1, cex=0.8)




#dev.off()


