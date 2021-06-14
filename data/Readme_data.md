# **Archivos:**
## **'full'** (Todos tienen las mismas propiedades)

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

**Nota**: El archivo voids tiene una columna mas:
- tipo: void tipo R (flag=0) y tipo S (flag=1)
*********************************************************************************************************
*********************************************************************************************************
## **'compact_for_vane'**
Es lo mismo Pero en la ultima fila se agrega el entorno (en character)

