 library(bootstrap); library(latex2exp); library(astrolibR); library(magicaxis)
#
#--------------------------------------------------------------
# Galaxias con line ade emision son activas: type spectral 1 ,2 ,3 (e(a),e(b) , e(c))
# El resto son galaxias pasivas: 4, 5, 6 ( k, k+a, a+k)
SS<-read.table("../data/compact_in_gg_m3_full")      ;S<-read.table("../data/compact_in_node_m3_full")
Gf<-read.table("../data/compact_in_field_m3_full");  FF<-read.table("../data/compact_in_filaments_m3_full")
HH<-read.table("../data/tab_gal_gru.dat"); VG<-read.table("../data/compact_in_voids_m3_full")
naranja <- rgb(1, 101/255, 0,1)

bin2=3

colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(HH)[c(1,2,3,4,5,6,7)] <-c("GId", "Nm", "RA", "Dec", "Redshift", "mag_r", "mag_g")

tcr<-c(S$tcr,SS$tcr,FF$tcr,Gf$tcr,VG$tcr)
sigma<-c(S$sigv,SS$sigv,FF$sigv,Gf$sigv,VG$sigv)
dij<-c(S$dij,SS$dij,FF$dij,Gf$dij,VG$dij)


mu_no<-S$tcr;  R_no<-S$rabs1 ; N_no=length(R_no)

mu_gg<-SS$tcr; R_gg<-SS$rabs1 ; N_gg=length(R_gg);

mu_fi<-FF$tcr; R_fi<-FF$rabs1 ; N_fi=length(R_fi)

mu_cp<- Gf$tcr ; R_cp<- Gf$rabs1 ; N_cp=length(R_cp)

mu_vv<-subset(VG$tcr, VG$V15==1) ; R_vv<-subset(VG$rabs1, VG$V15==1) ; N_vv=length(R_vv)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#	seleccion de galaxias miembros
for(i in 1:N_fi){
fil_gal_r<-subset(HH$mag_r, FF$igru[i]==HH$Gid )
}


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
  min_f=min(R_fi); max_f=max(R_fi)+0.1;    dM=(max_f-min_f)/(bin2)
  error_f<-vector("logical",bin2) ;   prop_f<-vector('logical',bin2)
  error2_f<-vector("logical",bin2) ;  count_f<-vector('logical',bin2)
  colmed_f<-vector("logical",bin2);       FL_f<-vector("logical",bin2)
  a_f<-vector("logical",N_fi)
  prop_f[1:bin2]=0

  for(i in 1:N_fi){
   ibin=as.integer((R_fi[i]-min_f)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    count_f[k]=count_f[k]+1
    prop_f[k]=prop_f[k]+mu_fi[i]
    a_f[i]<- mu_fi[i]
   }
  
  for(k in 1:bin2){
   colmed_f[k]=min_f+(dM/2.)*(2*k-1)
   FL_f[k]= (prop_f[k]/count_f[k])
   error2_f[k]=dM
  }


  boot=250; ff=0; jj=0
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
    a_gg[i]<- mu_gg[i]
   }




  boot=250; ff=0; jj=0
  count_g2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_gg[1:count_gg[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){mean(r)}  ;  results<-bootstrap(r,boot,theta)
#      for(i in 1:boot){
#      for(hh in 1:count_gg[k]){
#      jj=jj+results$thetastar[hh,i]
#      }
#      count_g2[i]<-jj/count_gg[k]
#     ff=0; jj=0
#     }
     error_gg[k]<-sd(results$thetastar) #sd(count_g2)
   }
#  for(k in 1:bin2){
#  if(error_gg[k] == 0){error_gg[k] = 0.02}
#  }

  for(k in 1:bin2){
   colmed_gg[k]=min_gg+(dM/2.)*(2*k-1)
   FL_gg[k]= (prop_gg[k]/count_gg[k])
#   error2_gg[k]=dM
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
    a_no[i]<- mu_no[i]
   }

  for(k in 1:bin2){
   colmed_no[k]=min_no+(dM/2.)*(2*k-1)+0.2
   FL_no[k]= (prop_no[k]/count_no[k])
   error2_no[k]=dM
  }



  boot=5000; ff=0; jj=0
  count_n2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_no[1:count_no[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){mean(r)}  ;  results<-bootstrap(r,boot,theta)
#     for(i in 1:boot){
#      for(hh in 1:count_no[k]){
#      jj=jj+results$thetastar[hh,i]
#      }
#      count_n2[i]<-jj/count_no[k]
#     ff=0; jj=0
#     }
     error_no[k]<-sd(results$thetastar)
   }
  for(k in 1:bin2){
 if(error_no[k] == 0){error_no[k] = 0.02}
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
    a_cp[i]<- mu_cp[i]
   }

  for(k in 1:bin2){
   colmed_cp[k]=min_cp+(dM/2.)*(2*k-1)
   FL_cp[k]= (prop_cp[k]/count_cp[k])
   error2_cp[k]=dM
  }



  boot=250; ff=0; jj=0
  count_n2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_cp[1:count_cp[k]] #;  rr<-aa_f[1:new_count_f[k],k]
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
     for(i in 1:boot){
      for(hh in 1:count_cp[k]){
      jj=jj+results$thetastar[hh,i]
      }
      count_n2[i]<-jj/count_cp[k]
     ff=0; jj=0
     }
     error_cp[k]<-sd(count_n2)
   }
  for(k in 1:bin2){
 if(error_cp[k] == 0){error_cp[k] = 0.02}
  }


#------------------------------------------------------------------------------------------------------------

  min_vv=min(R_vv); max_vv=max(R_vv);    dM=(max_vv-min_vv)/(bin2)
  error_vv<-vector("logical",bin2) ;   prop_vv<-vector('logical',bin2)
  error2_vv<-vector("logical",bin2) ;  
  colmed_vv<-vector("logical",bin2);       FL_vv<-vector("logical",bin2)
  a_vv<-vector("logical",N_vv)   
  count_vv<-vector(mode="logical",bin2) 
  count_vv[1:bin2]=0; prop_vv[1:bin2]=0

  for(i in 1:N_vv){
   ibin=as.integer((R_vv[i]-min_vv)/dM)+1
    if(ibin<1)   {ibin=1}
    if(ibin>bin2){ibin=bin2}
    k=ibin
    count_vv[k]=count_vv[k]+1
    prop_vv[k]=prop_vv[k]+mu_vv[i]
    a_vv[i]<- mu_vv[i]
   }

  for(k in 1:bin2){
   colmed_vv[k]=min_vv+(dM/2.)*(2*k-1)
   FL_vv[k]= (prop_vv[k]/count_vv[k])
   error2_vv[k]=dM
  }

  boot=2500; ff=0; jj=0
  count_c2<-vector(mode="logical",boot) 
  for(k in 1:bin2){
    r<-a_vv[1:count_vv[k]] 
    theta <- function(r){r}  ;  results<-bootstrap(r,boot,theta)
     for(i in 1:boot){
      for(hh in 1:count_vv[k]){
      jj=jj+results$thetastar[hh,i]
      }
      count_c2[i]<-jj/count_vv[k]
     ff=0; jj=0
     }
     error_vv[k]<-sd(count_c2)
   }
  for(k in 1:bin2){
 if(error_vv[k] == 0){error_vv[k] = 0.02}
  }


#print(count_no);print(count_gg); print(count_f);print(count_cp)
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
laby=TeX('$\\sigma$') 
labx=TeX('$M_r$')

par(family="serif")
par(cex.lab=0.7)       #Tamaño labels
par(cex.axis=0.6)      #Tamaño números en ejes.
par(lwd=2)
par(cex=2)
par(mgp=c(1.,0.1,0))      	        #margen labels
par(mar=c(4,4,4,4))


#plot(colmed_f,FL_f,type='l',ylim=c(0,0.15),xlim=c(-20,-23.),xlab=labx,ylab=TeX('$\\r_p$ \\[Kp h$^{-1}$ \\]'),
#col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n')
plot(colmed_f,FL_f,type='l',ylim=c(0,0.1),xlim=c(-20,-23.),xlab=labx,ylab=TeX('$\\t_cr$'),
col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n')
points(colmed_f,FL_f,pch=17,type='p',col=naranja,main='',lwd=3,asp=-5,xaxt='n',yaxt='n',cex=1)
polygon(c(colmed_f,rev(colmed_f)),c(FL_f+error_f*0.6,rev(FL_f-error_f*0.6)),col=rgb(1, 101/255, 0,0.5),border=NA)
#arrows(colmed_f,FL_f-error_f, colmed_f,FL_f+error_f,col=naranja,angle=90, code=3,length=0.1)

points(colmed_gg,FL_gg,type='l',col='magenta',main='',lwd=3)
points(colmed_gg,FL_gg,pch=18,type='p',col='magenta',main='')
polygon(c(colmed_gg,rev(colmed_gg)),c(FL_gg+error_gg,rev(FL_gg-error_gg)),col=rgb(1,153/255,1,0.5),border=NA)
#arrows(colmed_gg,FL_gg-error_gg, colmed_gg,FL_gg+error_gg,col='magenta',angle=90, code=3,length=0.1)

 points(colmed_no,FL_no,pch=16,type='l',col='red',main='',lwd=3)
 points(colmed_no,FL_no,pch=16,type='p',col='red',main='')
polygon(c(colmed_no,rev(colmed_no)),c(FL_no+error_no*0.6,rev(FL_no-error_no*0.6)),col=rgb(1,0,0,0.5),border=NA)
#arrows(colmed_no,FL_no-error_no, colmed_no,FL_no+error_no, col = 'red', angle=90, code=3,length=0.1)

 points(colmed_cp,FL_cp,type='l',col='darkblue',main='',lwd=3)
 points(colmed_cp,FL_cp,pch=20,type='p',col='darkblue',main='')
polygon(c(colmed_cp,rev(colmed_cp)),c(FL_cp+error_cp*0.6,rev(FL_cp-error_cp*0.6)),col=rgb(0,0,1,0.5),border=NA)
#arrows(colmed_cp,FL_cp-error_cp,colmed_cp,FL_cp+error_cp,col='darkblue',angle=90,code=3,length=0.1)


 points(colmed_vv,FL_vv,type='l',col='black',main='',lwd=3)
 points(colmed_vv,FL_vv,pch=21,type='p',col='black',main='')
 #arrows(colmed_vv,FL_vv-error_vv,colmed_vv,FL_vv+error_vv,col='black',angle=90,code=3,length=0.1)
polygon(c(colmed_vv,rev(colmed_vv)),c(FL_vv+error_vv,rev(FL_vv-error_vv)),col=rgb(0,0,0,0.4),border=NA)


legend(-20.,134,c("Node", "Filaments",'Groups','Field','Voids'),bty="n",lty=c(1,1,1,1,1), col=c('red',naranja,"magenta",'darkblue','black'),horiz=FALSE,inset=0,cex=0.5,pch=c(16,18,17,20,21))

 magaxis(1,majorn=5, minorn=5, tcl=0.3, ratio=0.5)
 magaxis(2,majorn=5, minorn=5, tcl=0.3, ratio=0.5)
 magaxis(3,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
 magaxis(4,majorn=5, minorn=5, tcl=0.3, ratio=0.5,labels=FALSE)
