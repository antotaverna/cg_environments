 library(bootstrap); library(latex2exp); library(astrolibR)
 library(magicaxis)
#--------------------------------------------------------------

#--------------------------------------------------------------
par(family="serif")
par(cex.lab=1.4)       #Tamaño labels
par(cex.main=1.1)        #Tamaño título
par(cex.axis=1.5)      #Tamaño números en ejes.
par(lwd=1)
par(cex=1)
par(mgp=c(1.8,0.4,0))      	        #margen labels
par(mar=c(4,4,4,4))
#----------------------------------------------------------------------------



SS<-read.table("compact_in_gg_m3_full")      ;S<-read.table("compact_in_node_m3_full")
Gf<-read.table("compact_in_field_m3_full");  FF<-read.table("compact_in_filaments_m3_full")
colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')
colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')
colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')
colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')


bin2=3
d="rabs1"
R_no<-S$rabs1  ;mu_no<-S$sigv  ; N_no=length(R_no)
R_gg<-SS$rabs1 ;mu_gg<-SS$sigv ; N_gg=length(R_gg)
R_fi<-FF$rabs1 ;mu_fi<-FF$sigv ; N_fi=length(R_fi)
R_cp<-Gf$rabs1 ;mu_cp<-Gf$sigv ; N_cp=length(R_cp)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
  min_f=min(R_fi); max_f=max(R_fi);    dM=(max_f-min_f)/(bin2)
  error_f<-vector("logical",bin2) ;   prop_f<-vector('logical',bin2)
  error2_f<-vector("logical",bin2);   count_f<-vector('logical',bin2)
  colmed_f<-vector("logical",bin2);       FL_f<-vector("logical",bin2)
  a_f<-vector("logical",N_fi)
  prop_f[1:bin2]=0

  for(i in 1:N_fi){
   ibin=as.integer((R_fi[i]-min_f)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
#    if (R_fi[i] <= -18.0 & R_fi[i] >=-21.4){k=1}
#    if (R_fi[i] <= -21.4 & R_fi[i] >=-21.8){k=2}
#    if (R_fi[i] <= -21.8 & R_fi[i] >=-24.0){k=3}
    count_f[k]=count_f[k]+1
    prop_f[k]=prop_f[k]+mu_fi[i]
    a_f[count_f[k]]<- mu_fi[i]
   }

  for(k in 1:bin2){
   colmed_f[k]=min_f+(dM/2.)*(2*k-1)
   FL_f[k]= (prop_f[k]/count_f[k])
   error2_f[k]=dM
  }


  boot=200; ff=0; jj=0
  count_ff<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_f[1:count_f[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
#    theta <- function(rr){rr};results2<-bootstrap(rr,boot,theta)
     for(i in 1:boot){
      for(hh in 1:count_f[k]){
      jj=jj+results$thetastar[hh,i]
      }
      count_ff[i]<-jj/count_f[k]
     ff=0; jj=0
     }
     error_f[k]<-sd(count_ff)
   }
  for(k in 1:bin2){
 if(error_f[k] == 0){error_f[k] = 0.02}
  }

#------------------------------------------------------------------------------------------------------------

  min_gg=min(R_gg); max_gg=max(R_gg);    dM=(max_gg-min_gg)/(bin2)
  error_gg<-vector("logical",bin2) ;   prop_gg<-vector('logical',bin2)
  error2_gg<-vector("logical",bin2) ;  count_g2<-vector('logical',bin2)
  colmed_gg<-vector("logical",bin2);       FL_gg<-vector("logical",bin2)
  a_gg<-vector("logical",N_gg)   
  count_gg<-vector(mode="logical",bin2) 
  count_gg[1:bin2]=0; prop_gg[1:bin2]=0

  for(i in 1:N_gg){
   ibin=as.integer((R_gg[i]-min_gg)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    	count_gg[k]=count_gg[k]+1
    prop_gg[k]=prop_gg[k]+mu_gg[i]
    a_gg[count_gg[k]]<- mu_gg[i]
   }
  for(k in 1:bin2){
   colmed_gg[k]=min_gg+(dM/2.)*(2*k-1)
   FL_gg[k]= (prop_gg[k]/count_gg[k])
   error2_gg[k]=dM
  }



  boot=250; ff=0; jj=0
  count_g2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_gg[1:count_gg[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
     for(i in 1:boot){
      for(hh in 1:count_gg[k]){
      jj=jj+results$thetastar[hh,i]
      }
      count_g2[i]<-jj/count_gg[k]
     ff=0; jj=0
     }
     error_gg[k]<-sd(count_g2)
   }
  for(k in 1:bin2){
 if(error_gg[k] == 0){error_gg[k] = 0.02}
  }

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
  min_no=min(R_no); max_no=max(R_no);    dM=(max_no-min_no)/(bin2)
  error_no<-vector("logical",bin2) ;   prop_no<-vector('logical',bin2)
  error2_no<-vector("logical",bin2) ;  
  colmed_no<-vector("logical",bin2);       FL_no<-vector("logical",bin2)
  a_no<-vector("logical",N_no)   
  count_no<-vector(mode="logical",bin2) 
  count_no[1:bin2]=0; prop_no[1:bin2]=0

  for(i in 1:N_no){
   ibin=as.integer((R_no[i]-min_no)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    count_no[k]=count_no[k]+1
    prop_no[k]=prop_no[k]+mu_no[i]
    a_no[count_no[k]]<- mu_no[i]
   }

  for(k in 1:bin2){
   colmed_no[k]=min_no+(dM/2.)*(2*k-1)
   FL_no[k]= (prop_no[k]/count_no[k])
   error2_no[k]=dM
  }



  boot=250; ff=0; jj=0
  count_n2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_no[1:count_no[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
     for(i in 1:boot){
      for(hh in 1:count_no[k]){
      jj=jj+results$thetastar[hh,i]
      }
      count_n2[i]<-jj/count_no[k]
     ff=0; jj=0
     }
     error_no[k]<-sd(count_n2)
   }
  for(k in 1:bin2){
 if(error_no[k] == 0){error_no[k] = 10}
  }


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
  min_cp=min(R_cp); max_cp=max(R_cp);    dM=(max_cp-min_cp)/(bin2)
  error_cp<-vector("logical",bin2) ;   prop_cp<-vector('logical',bin2)
  error2_cp<-vector("logical",bin2) ;  
  colmed_cp<-vector("logical",bin2);       FL_cp<-vector("logical",bin2)
  a_cp<-vector("logical",N_cp)   
  count_cp<-vector(mode="logical",bin2) 
  count_cp[1:bin2]=0; prop_cp[1:bin2]=0

  for(i in 1:N_cp){
   ibin=as.integer((R_cp[i]-min_cp)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    count_cp[k]=count_cp[k]+1
    prop_cp[k]=prop_cp[k]+mu_cp[i]
    a_cp[count_cp[k]]<- mu_cp[i]
   }

  for(k in 1:bin2){
   colmed_cp[k]=min_cp+(dM/2.)*(2*k-1)
   FL_cp[k]= (prop_cp[k]/count_cp[k])
   error2_cp[k]=dM
  }

  boot=250; ff=0; jj=0
  count_c2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_cp[1:count_cp[k]] 
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
     for(i in 1:boot){
      for(hh in 1:count_cp[k]){
      jj=jj+results$thetastar[hh,i]
      }
      count_c2[i]<-jj/count_cp[k]
     ff=0; jj=0
     }
     error_cp[k]<-sd(count_c2)
   }
  for(k in 1:bin2){
 if(error_cp[k] == 0){error_cp[k] = 10}
  }


print(count_no);print(count_gg); print(count_f);print(count_cp)
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
laby=TeX('$\\sigma$') 
labx=TeX('$M_r$')


xm <- colmed_f ;ym <- FL_f ; sd <- error_f
plot(colmed_f,FL_f,type='l',ylim=c(0,600),xlim=c(-24,-19),xlab=labx,ylab=laby,
col='darkgoldenrod1',main='',lwd=3,asp=-5,xaxt='n',yaxt='n')
points(colmed_f,FL_f,pch=17,type='p',col='darkgoldenrod1',main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(xm,rev(xm)),c(ym+sd,rev(ym-sd)),col=rgb(1, 101/255, 0,0.5),border=NA)


xm <- colmed_gg ;ym <- FL_gg ; sd <- error_gg
points(colmed_gg,FL_gg,type='l',col='magenta',main='',lwd=3)
points(colmed_gg,FL_gg,pch=18,type='p',col='magenta',main='')
polygon(c(xm,rev(xm)),c(ym+sd,rev(ym-sd)),col=rgb(1, 0.2, 0.8,0.5),border=NA) 

xm <- colmed_no ;ym <- FL_no ; sd <- error_no
points(colmed_no,FL_no,pch=16,type='l',col='red',main='',lwd=3)
points(colmed_no,FL_no,pch=16,type='p',col='red',main='')
polygon(c(xm,rev(xm)),c(ym+sd,rev(ym-sd)),col=rgb(1,0,0,0.5),border=NA) 


xm <- colmed_cp ;ym <- FL_cp ; sd <- error_cp
points(colmed_cp,FL_cp,pch=16,type='l',col='darkblue',main='',lwd=3)
points(colmed_cp,FL_cp,pch=20,type='p',col='darkblue',main='')
polygon(c(xm,rev(xm)),c(ym+sd,rev(ym-sd)),col= rgb(0,0,1,0.5),border=NA)


legend('topleft',c("Node",'Groups', "Filaments",'Field'),bty="n",lty=c(1,1,1,1), col=c('red',"magenta","darkgoldenrod1",'darkblue'),horiz=FALSE,inset=0,cex=1.2,pch=c(16,17,18,20))

 magaxis(1,majorn=5, minorn=5, tcl=0.5, ratio=0.5)
 magaxis(2,majorn=5, minorn=5, tcl=0.5, ratio=0.5)
 magaxis(3,majorn=5, minorn=5, tcl=0.5, ratio=0.5,labels=FALSE)
 magaxis(4,majorn=5, minorn=5, tcl=0.5, ratio=0.5,labels=FALSE)






