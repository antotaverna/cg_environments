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

*********************************************************************************************************
*********************************************************************************************************
# Plots
## box_plot.r
Realiza el box_plot de propiedaes de compactos (boxplot.pdf).  


## brillo_superficial.r
Realiza comparaciones de propiedades Y en funcion de X. En el overleafe hay dos ejemplos 
(sigma_vs_mr.pdf, r_p_vs_rRmag.pdf).


## gap_mag.r
Box plot de la diferencia de las dos galaxias mas brillante (dif_m2_m1.pdf).


## color_gal_in_cg.r
Trabaja con las galaxias en grupos compacto. Calcula la distribucion de color (selección de 
galaxias pasivas). Fracción de galaxias pasivas (colorgr_galcg.pdf; frac_pas_all_galaxies.pdf).


## all_frac_vs_mag.R
Trabaja con la fraccion de galaxias pasivas (La misma que en color_gal_in_cg.r)
                   (frac_red_gal_mag.pdf)


## vol_lim.r
Redshift vs Mabs. volume-limited sample (vol_lim.pdf).


## aitoff.R
aitoff (aitoff.pdf).
