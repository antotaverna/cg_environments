# **Archivos:**
## **'tab_gal_grupos.dat'**

- GId: Compact Group ID. 
- Nm: Galaxy index.
- RA: Right Ascension (J2000).
- Dec: Declination (J2000).
- Redshift: Galaxy redshift corrected to the CMB rest-frame.
- mag_r : r-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.
- mag_g : g-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.
- mag_u : u-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.
- galID : galaxy ID-SDSS.  
- r50: petrosian radius
- r90: petrosian radius
- abs_r : r-band k-corrected absolute magnitud.
- abs_g : g-band k-corrected absolute magnitud.
- abs_u : u-band k-corrected absolute magnitud.
- glx_pasiva/red : flag=1 (abs_u-abs_r.gt.C(abs_r)).
- glx_early-type : flag=1 (r90/r50.gt.2.5).

## **'tab_grupos.dat'**

- i_gru: ID del grupo compacto
- n_mi:  Numeros de miembros 
- Ra:    Asención recta
- dec:   Declinación
- z:     Redshift medio del CG
- radio_mins: Tita_g /2 
- mu:    Brillo superficial medio (banda r).
- sigv:  Dispersión de velocidad
- rmag_bri: Magnitud aparente de la galaxia mas brillate  
- flag: Entero que señala que un Cg este o no contaminado (potencialmente). 
- d_ij: Distancia mediana entre miembros. 
- tcr:  Tiemto de crossing.
- rp:   Radio proyectado (radio mins en Kpc).
- rabs1: Magnitud absoluta de la galaxia mas brillante.
- rabs2: Magnitud absoluta de la 2da galaxia mas brillante. 


## **compact_in_xxx_full** (Todos tienen las mismas propiedades)

- 'compact_in_gg_m3_full'
- 'compact_in_node_m3_full'
- 'compact_in_field_m3_full'
- 'compact_in_filaments_m3_full'
- 'compact_in_voids_m3_full'

**Columnas**
- i_gru: ID del grupo compacto
- n_mi:  Numeros de miembros 
- Ra:    Asención recta
- dec:   Declinación
- z:     Redshift medio del CG
- radio_mins: Tita_g /2 
- mu:    Brillo superficial medio (banda r).
- sigv:  Dispersión de velocidad
- rmag_bri: Magnitud aparente de la galaxia mas brillate  
- flag: Entero que señala que un Cg este o no contaminado (potencialmente).
- rabs1: Magnitud absoluta de la galaxia mas brillante. 
- d_ij: Distancia mediana entre miembros. 
- tcr:  Tiemto de crossing.
- rp:   Radio proyectado (radio mins en Kpc).
- rabs2: Magnitud absoluta de la 2da galaxia mas brillante. 

**Nota**: El archivo voids tiene dos columnas mas:
- tipo: void tipo R (flag=0) y tipo S (flag=1)
- estado: void adentro (flag=-1) y void borde (flag=0)
*********************************************************************************************************
*********************************************************************************************************
## **'compact_in_all'**
Es lo mismo que el archivo compact_in_xxx Pero en la ultima fila se agrega el entorno (en character)

## **'compact_in_all_estado.csv'**
Es lo mismo que el archivo compact_in_all Pero en la ultima fila se cambian las etiquetas de void voidIN y voidED

## **tempel_limpio.dat**: DR12 SDSS sample of galaxies
- galID   :ID galaxia (SDSS)
- ra[rad] :Asención recta
- dec[rad]:Declinacion
- zCMB   :Galaxy redshift corrected to the CMB rest-frame.
- rextAB :r-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.
- gextAB :g-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.

## **tempel_Mr.dat**: 
- Se le agregan 3 columnas al tempel_limpio.dat
- rabs :r-band k-corrected absolute magnitude.
- gabs :g-band k-corrected absolute magnitude.
- Menv :absolute magntude "r", Menv= mlim - 25 -5*log(ld) -kr, (mlim=17.77) .

