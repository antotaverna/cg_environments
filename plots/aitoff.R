#source("aitoff.R")

CG <- read.table("/home/ataverna/Documentos/proyectos/cg_in_different_environments/catalogos/anto_gitano_vane/muestras_filtradas/cg/cgCV.dat",head=FALSE)

VR <- read.table('/home/ataverna/Documentos/proyectos/cg_in_different_environments/catalogos/anto_gitano_vane/muestras_filtradas/void/voids_filtrados_r.dat',head=FALSE)

VS <- read.table('/home/ataverna/Documentos/proyectos/cg_in_different_environments/catalogos/anto_gitano_vane/muestras_filtradas/void/voids_filtrados_s.dat',head=FALSE)

GXS_V <- read.table('/home/ataverna/Documentos/proyectos/cg_in_different_environments/catalogos/DR12_tempel.dat',head=FALSE)
#read.table('gxs_v.dat',head=FALSE) -> GXS_V


rac_cg <- CG$V3
dec_cg <- CG$V4
rac_vr <- VR$V2
dec_vr <- VR$V3
rac_vs <- VS$V2
dec_vs <- VS$V3





###############################
library(dplyr)

#rac_gxs <- GXS_V$V1
#dec_gxs <- GXS_V$V2
rac_gxs_all <- GXS_V$V2*180/pi
dec_gxs_all <- GXS_V$V3*180/pi

df=data.frame(rac_gxs_all,dec_gxs_all)
#sample_n(df,4)
new_df=sample_frac(df,0.05)

dim(new_df); nrow(new_df); ncol(new_df)

rac_gxs <- new_df[,1]
dec_gxs <- new_df[,2]


######### PLOT #################
library(SPADAR)
#svg("aitoff_R.svg")
#pdf("aitoff_R.pdf")



createAllSkyGridChart(longitude = c(0, 45, 90, 135, 180, 225, 270, 315, 360),
                      latitude = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75), 
		      mainGrid = "equatorial",eqCol = "gray0", eclCol = "blue", 
		      galCol = "green", eqLty = 1, eclLty = 2, galLty = 3,
                      eqLwd = 1, eclLwd = 1, galLwd = 1, eqDraw = TRUE, eclDraw = FALSE,
		      galDraq = FALSE, projname = "aitoff", projparam = NULL, 
		      projorient = NULL, npoints = 100, overplot = FALSE, addLab = TRUE, 
		      label.cex = 0.6, main = NULL)
#galaxias tempel
overplotScatterPlotInAllSkyGridChart(rac_gxs, dec_gxs, pointcol = 'grey', # azure4
				     dataCoordSys = "equatorial",
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL, 
				     pch='.',pointsize =0.02)
#dorados CG
overplotScatterPlotInAllSkyGridChart(rac_cg, dec_cg, pointcol = "red", #darkgoldenrod1,
                                     dataCoordSys = "equatorial", 
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL,
				     pch=16,pointsize =0.5)
#verdes Voids tipo R
overplotScatterPlotInAllSkyGridChart(rac_vr, dec_vr, pointcol = 'darkgreen',
                                     dataCoordSys = "equatorial",
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL,
				     pch=1,pointsize =0.5)
				     #pch='.',pointsize =4)
#azul Voids tipo S
overplotScatterPlotInAllSkyGridChart(rac_vs, dec_vs, pointcol = 'blue',
                                     dataCoordSys = "equatorial", 
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL, 
				     pch=1,pointsize =0.5)
				     #pch='.',pointsize =4)


#dev.off()


#########################################################
#########################################################
########### aitoff CG in XX #############################
#########################################################
#########################################################

SS<-read.table("../muestras_finales/compact_in_gg_m3_full")
S<-read.table("../muestras_finales/compact_in_node_m3_full")
Gf<-read.table("../muestras_finales/compact_in_field_m3_full")
FF<-read.table("../muestras_finales/compact_in_filaments_m3_full")
VG<-read.table("../muestras_finales/compact_in_voids_m3_full")

gal<-read.table("../catalogos/DR12_tempel.dat")


colnames(SS)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('igru','nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(Gf)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(S)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(FF)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')

colnames(VG)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c("igru",'nmi','alcm','delcm','zmedian','radio_mins','mu','sigv','rmag_bri','flag','rabs1','dij','tcr','rp')


al_no<-S$alcm;  del_no<-S$delcm

al_gg<-SS$alcm; del_gg<-SS$delcm

al_fi<-FF$alcm; del_fi<-FF$delcm

al_cp<-Gf$alcm; del_cp<-Gf$delcm

al_vv<-VG$alcm; del_vv<-VG$delcm

rac_gxs_all <- gal$V2*180/pi
dec_gxs_all <- gal$V3*180/pi

df=data.frame(rac_gxs_all,dec_gxs_all)
#sample_n(df,4)
new_df=sample_frac(df,0.05)

dim(new_df); nrow(new_df); ncol(new_df)

rac_gxs <- new_df[,1]
dec_gxs <- new_df[,2]



######### PLOT #################
library(SPADAR)
#svg("aitoff_R.svg")
#pdf("aitoff_R.pdf")



createAllSkyGridChart(longitude = c(0, 45, 90, 135, 180, 225, 270, 315, 360),
                      latitude = c(-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75), 
		      mainGrid = "equatorial",eqCol = "gray0", eclCol = "blue", 
		      galCol = "green", eqLty = 1, eclLty = 2, galLty = 3,
                      eqLwd = 1, eclLwd = 1, galLwd = 1, eqDraw = TRUE, eclDraw = FALSE,
		      galDraq = FALSE, projname = "aitoff", projparam = NULL, 
		      projorient = NULL, npoints = 100, overplot = FALSE, addLab = TRUE, 
		      label.cex = 0.6, main = NULL)
#galaxias tempel
overplotScatterPlotInAllSkyGridChart(rac_gxs, dec_gxs, pointcol = 'grey', # azure4
				     dataCoordSys = "equatorial",
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL, 
				     pch='.',pointsize =0.02)

# CG in NODES
overplotScatterPlotInAllSkyGridChart(al_no, del_no, pointcol = "red", #darkgoldenrod1,
                                     dataCoordSys = "equatorial", 
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL,
				     pch=16,pointsize =0.7)

# CG in LOOSE
overplotScatterPlotInAllSkyGridChart(al_gg, del_gg, pointcol = 'magenta',
                                     dataCoordSys = "equatorial",
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL,
				     pch=0,pointsize =0.5)
				     #pch='.',pointsize =4)

# CG in FILAMENT
overplotScatterPlotInAllSkyGridChart(al_fi, del_fi, pointcol = 'orange',
                                     dataCoordSys = "equatorial", 
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL, 
				     pch=5,pointsize =0.5)
				     #pch='.',pointsize =4)

# CG in FIELD
overplotScatterPlotInAllSkyGridChart(al_cp, del_cp, pointcol = 'blue',
                                     dataCoordSys = "equatorial", 
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL, 
				     pch=4,pointsize =0.5)
				     #pch='.',pointsize =4)

# CG in VOIDS
overplotScatterPlotInAllSkyGridChart(al_vv, del_vv, pointcol = 'black',
                                     dataCoordSys = "equatorial",
				     mainGrid = "equatorial",
                                     projname = "aitoff", projparam = NULL,
				     pch=11,pointsize =0.5)
				     #pch='.',pointsize =4)



#dev.off()


