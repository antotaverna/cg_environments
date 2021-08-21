# source("GMM_mclust.r")
library(MASS)
library("mclust")
library(ggplot2)
require(mixtools)
library(AdaptGauss)

   
#############################################
#Muestra galaxias tempel
#############################################
glx<-read.table("../data/tempel_Mr.dat",sep="")
colnames(glx)[c(1,2,3,4,5,6,7,8)] <- c("id","al","del","z","mr","mg","rabs","gabs")

r=subset(glx$rabs,glx$mg>0 & glx$gabs>-26)
g=subset(glx$gabs,glx$mg>0 & glx$gabs>-26)

color <- subset(g - r, g-r>0 & g-r<1.1)
rr <- subset(r, g-r>0 & g-r<1.1)
anyNA(color)

#Mr=[-23,-22)
color1= subset(color, rr>=-23. & rr< -22.)
#Mr=[-22,-21.5)
color2= subset(color,rr>=-22. & rr< -21.5)
#Mr=[-21.5,-21)
color3= subset(color,rr>=-21.5 & rr< -21.)
#Mr=[-21.,-20.5)
color4= subset(color,rr>=-21. & rr< -20.5)
#Mr=[-20.5.,-20.)
color5= subset(color,rr>=-20.5 & rr< -20.)
#Mr=[-20.,-19.5)
color6= subset(color,rr>=-20. & rr< -19.5)
#Mr=[-19.5.,-19.)
color7= subset(color,rr>=-19.5 & rr< -19.)
#Mr=[-19..,-18.)
color8= subset(color,rr>=-19. & rr< -18.)

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


#summary(gmm, parameters = TRUE)
#gmm$classification --> a que gassiana pertenece
############################################
#    PLOT 
############################################

#--------------------------------------
plot_col <- function(xx,gaus,a){
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
legend("topleft", legend=c(a), cex=0.8)
dd=data.frame(a,ll)
#write.table(dd,"gmm_intersect.dat",sep=" ",row.names=F,col.names=F, append = TRUE)
xtick<-seq(0, 10, by=5)
axis(side=1, at=c(0.,0.2,0.4,0.6,0.8,1.0,1.2), labels = FALSE)
}
#------------------------------------

#cairo_pdf("color_tempel_gmm.pdf")
#cairo_ps("color_tempel_gmm.eps")
par(mfrow=c(3,3))
par(mar = c(2, 2, 2, 2))
par(oma = c(5, 5, 2, 2))
#xx=color1
#gaus=gmm1

plot_col(color,gmm,"[-23,-17]")
plot_col(color1,gmm1,"[-23,-22]")
plot_col(color2,gmm2,"[-22,-21.5]")
plot_col(color3,gmm3,"[-21.5,-21]")
plot_col(color4,gmm4,"[-21,-20.5]")
plot_col(color5,gmm5,"[-20.5,-20]")
plot_col(color6,gmm6,"[-20,-19.5]")
plot_col(color7,gmm7,"[-19.5,-19]")
plot_col(color8,gmm8,"[-19,-18]")

mtext(expression(paste(g,'\u2212',r," ","[mag]")), side = 1, cex = 1, line = 3.2, col = "black", outer = TRUE)
mtext(expression(paste("Density")), side = 2, cex = 1, line = 3.2, col = "black", outer = TRUE)


#dev.off()

#par(ask=TRUE)




#---------------------------------
#cairo_ps("gmm.eps")
#cairo_pdf("gmm.pdf")
#dev.off()





