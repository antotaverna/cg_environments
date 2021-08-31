# source("GMM_mclust.r")
library(MASS)
library("mclust")
library(ggplot2)
#require(mixtools)
library(AdaptGauss)

   
#############################################
#Muestra galaxias tempel
#############################################
g<-read.table("../../tempel_Mr.dat",sep="")
#"id","al","del","z","mr","mg","rabs","gabs"

rma=subset(g$V7, g$V8 > -27. & g$V8 < -15.)
gma=subset(g$V8, g$V8 > -27. & g$V8 < -15.)
c = gma - rma

color <- subset(c, c > 0. & c < 1.1)
rr <- subset(rma, c > 0. & c < 1.1)
anyNA(color)
#bin=seq(-22.5,-17,0.35)
bin=seq(-22.5,-18.,0.3)


color1= subset(color, rr>=bin[1] & rr< bin[2])
color2= subset(color, rr>=bin[2] & rr< bin[3])
color3= subset(color, rr>=bin[3] & rr< bin[4])
color4= subset(color, rr>=bin[4] & rr< bin[5])
color5= subset(color, rr>=bin[5] & rr< bin[6])
color6= subset(color, rr>=bin[6] & rr< bin[7])
color7= subset(color, rr>=bin[7] & rr< bin[8])
color8= subset(color, rr>=bin[8] & rr< bin[9])
color9= subset(color, rr>=bin[9] & rr< bin[10])
color10= subset(color, rr>=bin[10] & rr< bin[11])
color11= subset(color, rr>=bin[11] & rr< bin[12])
color12= subset(color, rr>=bin[12] & rr< bin[13])
color13= subset(color, rr>=bin[13] & rr< bin[14])
color14= subset(color, rr>=bin[14] & rr< bin[15])
color15= subset(color, rr>=bin[15] & rr< bin[16])

############################################
#  2 Gaussians 
############################################
#gmm1 <- densityMclust(color,G=2,modelNames='VVV',initialization = set.seed(0))
gmm <- densityMclust(color,G=2,initialization = set.seed(0))
gmm1 <- densityMclust(color1,G=2,initialization = set.seed(0))
gmm2 <- densityMclust(color2,G=2,initialization = set.seed(0))
gmm3 <- densityMclust(color3,G=2,initialization = set.seed(0))
gmm4 <- densityMclust(color4,G=2,initialization = set.seed(0))
gmm5 <- densityMclust(color5,G=2,initialization = set.seed(0))
gmm6 <- densityMclust(color6,G=2,initialization = set.seed(0))
gmm7 <- densityMclust(color7,G=2,initialization = set.seed(0))
gmm8 <- densityMclust(color8,G=2,initialization = set.seed(0))
gmm9 <- densityMclust(color9,G=2,initialization = set.seed(0))
gmm10 <- densityMclust(color10,G=2,initialization = set.seed(0))
gmm11 <- densityMclust(color11,G=2,initialization = set.seed(0))
gmm12 <- densityMclust(color12,G=2,initialization = set.seed(0))
gmm13 <- densityMclust(color13,G=2,initialization = set.seed(0))
gmm14 <- densityMclust(color14,G=2,initialization = set.seed(0))
gmm15 <- densityMclust(color15,G=2,initialization = set.seed(0))


#summary(gmm, parameters = TRUE)
#gmm$classification --> a que gassiana pertenece
############################################
#    PLOT 
############################################

#--------------------------------------
plot_col <- function(xx,gaus,a,bin1,bin2){
br <- seq(min(xx), max(xx), length = 100)
#plot(gaus, what = "density", data = xx, breaks = br,xlab="",xaxt='n',yaxt='n')
plot(gaus, what = "density",  breaks = br,xlab="",lty=1)
h=hist(xx,plot=FALSE)
points(h$mids,h$density,type="s", col = "black")

#if(mark==0){plot(gaus, what = "density", breaks = br,xlab="",xaxt='n')}
#if(mark==1){plot(gaus, what = "density", breaks = br,xlab="",yaxt='n')}
#if(mark==2){plot(gaus, what = "density", breaks = br,xlab="",xaxt='n',yaxt='n')}
#if(mark==3){plot(gaus, what = "density", breaks = br,xlab="")}

x <- seq(min(xx)-diff(range(xx))/10, max(xx)+diff(range(xx))/10, length = 200)
cdens <- predict(gaus, x, what = "cdens")
cdens <- t(apply(cdens, 1, function(d) d*gaus$parameters$pro))
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 2:2, col = c("red","blue"),xaxt='n',yaxt='n')

m1=gaus$parameters$mean[1]; sd1=sqrt(gaus$parameters$variance$sigmasq[1]); p1=gaus$parameters$pro[1]
m2=gaus$parameters$mean[2]; sd2=sqrt(gaus$parameters$variance$sigmasq[2]); p2=gaus$parameters$pro[2]

#inteseccion de dos gaussianas
intersec=Intersect2Mixtures(m1,sd1,p1,m2,sd2,p2)
ll=intersec$CutX
abline(v=ll,col='magenta')
#legend("topleft", legend=c("[a,b]"), cex=0.8)
legend("topleft", legend=c(a), cex=0.7,bty="n")

lll=round(ll,digits=4)
dd=data.frame(bin1,bin2,lll)
#write.table(dd,"gmm_intersect.dat",sep="    ", dec = ".", quote = FALSE,row.names=F,col.names=F, append = TRUE)
xtick<-seq(0, 10, by=5)
axis(side=1, at=c(0.,0.2,0.4,0.6,0.8,1.0,1.2), labels = FALSE)
}
#------------------------------------

#cairo_pdf("color_tempel_gmm.pdf")
#cairo_ps("color_tempel_gmm.eps")
par(mfrow=c(4,4))
par(mar = c(2, 2, 2, 2))
par(oma = c(5, 5, 2, 2))
#xx=color1
#gaus=gmm1

plot_col(color,gmm,"[-23,-17]",-23,-17)
plot_col(color1,gmm1,as.character(paste("[",bin[1],",",bin[2],"]")),bin[1],bin[2])
plot_col(color2,gmm2,as.character(paste("[",bin[2],",",bin[3],"]")),bin[2],bin[3])
plot_col(color3,gmm3,as.character(paste("[",bin[3],",",bin[4],"]")),bin[3],bin[4])
plot_col(color4,gmm4,as.character(paste("[",bin[4],",",bin[5],"]")),bin[4],bin[5])
plot_col(color5,gmm5,as.character(paste("[",bin[5],",",bin[6],"]")),bin[5],bin[6])
plot_col(color6,gmm6,as.character(paste("[",bin[6],",",bin[7],"]")),bin[6],bin[7])
plot_col(color7,gmm7,as.character(paste("[",bin[7],",",bin[8],"]")),bin[7],bin[8])
plot_col(color8,gmm8,as.character(paste("[",bin[8],",",bin[9],"]")),bin[8],bin[9])
plot_col(color9,gmm9,as.character(paste("[",bin[9],",",bin[10],"]")),bin[9],bin[10])
plot_col(color10,gmm10,as.character(paste("[",bin[10],",",bin[11],"]")),bin[10],bin[11])
plot_col(color11,gmm11,as.character(paste("[",bin[11],",",bin[12],"]")),bin[11],bin[12])
plot_col(color12,gmm12,as.character(paste("[",bin[12],",",bin[13],"]")),bin[12],bin[13])
plot_col(color13,gmm13,as.character(paste("[",bin[13],",",bin[14],"]")),bin[13],bin[14])
plot_col(color14,gmm14,as.character(paste("[",bin[14],",",bin[15],"]")),bin[14],bin[15])
plot_col(color15,gmm15,as.character(paste("[",bin[15],",",bin[16],"]")),bin[15],bin[16])

mtext(expression(paste(g,'\u2212',r," ","[mag]")), side = 1, cex = 1, line = 3.2, col = "black", outer = TRUE)
mtext(expression(paste("Density")), side = 2, cex = 1, line = 3.2, col = "black", outer = TRUE)


#dev.off()





#---------------------------------
#cairo_ps("gmm.eps")
#cairo_pdf("gmm.pdf")
#dev.off()





