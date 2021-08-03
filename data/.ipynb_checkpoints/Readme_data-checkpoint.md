# **Archivos:**
## **'tab_gal_grupos.dat'**

- GId: Compact Group ID. 
- Nm: Galaxy index.
- RA: Right Ascension (J2000).
- Dec: Declination (J2000).
- Redshift: Galaxy redshift corrected to the CMB rest-frame.
- mag_r : r-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.
- mag_g : g-band observer-frame model apparent magnitud corrected for galactic extinction in AB system.
- abs_r : r-band k-corrected absolute magnitud
- abs_g : g-band k-corrected absolute magnitud
- glx_pasiva : flag=1 (abs_g-abs_r.gt.0.72 - 0.03 \times (abs_r+20))

## **'tab_grupos.dat' y xxx_full** (Todos tienen las mismas propiedades)

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

