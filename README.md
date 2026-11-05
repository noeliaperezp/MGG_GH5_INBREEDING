
# MGG GH5 - ANÁLISIS DE CONSANGUINIDAD

La consanguinidad (consecuencia del apareamiento entre individuos emparentados) es inevitable en poblaciones pequeñas en ausencia de migración, incluso cuando los apareamientos ocurren al azar. La consanguinidad suele ir asociada a una reducción de la eficacia media de las poblaciones, en comparación con poblaciones no consanguíneas, fenómeno conocido como depresión consanguínea. Esta reducción de la eficacia se debe principalmente al aumento de la homocigosis de mutaciones deletéreas (parcialmente) recesivas, que en poblaciones no consanguíneas permanecen enmascaradas en heterocigosis. La depresión consanguínea tiene por tanto importantes consecuencias para la viabilidad y resiliencia de las poblaciones.

La depresión consanguínea se puede estimar a partir de la pendiente de la regresión lineal de los valores fenotípicos de los individuos en función de sus coeficientes de consanguinidad. El grado de consanguinidad de un individuo se cuantifica mediante el coeficiente de consanguinidad *F* (Wright, 1922), que representa la probabilidad de que dos alelos en un locus sean idénticos por descendencia (IBD). Tradicionalmente, este coeficiente se ha estimado a partir de genealogías, lo que conlleva una serie de limitaciones: las genealogías suelen ser incompletas o inexistentes y solo reflejan expectativas teóricas, por lo que pueden diferir del nivel real de consanguinidad de cada individuo. Actualmente, la gran disponibilidad de marcadores moleculares (SNPs) permite estimar la consanguinidad de un individuo sin necesidad de datos genealógicos. A estas medidas se les conoce como consanguinidad molecular o genómica, y proporcionan estimaciones más precisas de la fracción del genoma que realmente es idéntica por descendencia, por lo que constituyen una herramienta muy útil para detectar depresión consanguínea. 

Existen distintos estimadores genómicos de *F*. Algunos se basan en el análisis SNP por SNP, que dependen de las frecuencias alélicas (e.g., el coeficiente *F* de Yang; Yang *et al.*, 2011), y otros se basan en la detección de Runs of Homozygosity (ROH), que son regiones del genoma que son homocigotas para la gran mayoría de las bases nucleotídicas (McQuillan *et al.*, 2008). Los métodos basados en el análisis SNP por SNP suelen utilizar las frecuencias de la generación actual, ya que, por lo general, se desconocen las frecuencias de generaciones anteriores (población de partida o población de referencia). En estos casos, estas estimaciones reflejan desviaciones de las frecuencias esperadas en el equilibrio de Hardy-Weinberg o correlaciones entre alelos, por lo que pueden adoptar valores negativos. Las estimaciones basadas en ROH, en cambio, proporcionan medidas de IBD. Aún existe debate sobre qué estimador refleja mejor la depresión consanguínea, pero tanto las medidas basadas en SNP como las basadas en ROH han demostrado ser eficaces para detectarla. En poblaciones donde no hay información histórica ni genealogías, las estimaciones genómicas de *F* ofrecen una herramienta práctica y potente para estudiar la consanguinidad y sus efectos sobre la eficacia.

En esta sesión práctica aprenderemos a calcular la consanguinidad genómica con diversos estimadores y a estimar la depresión consanguínea a partir de ellos. Para ello, dividiremos la sesión en tres secciones:
1. Estimación de la consanguinidad genómica en poblaciones simuladas:
   - Simulación de las poblaciones de estudio con el programa SLiM:
      - Población de gran tamaño con panmixia.
      - Población de gran tamaño que sufrió un cuello de botella hace 50 generaciones.
      - Población donde existe preferencia por apareamientos entre parientes.  
   - Cálculo de varios estimadores de consanguinidad con el programa *PLINK*
2. Estimación de la consanguinidad genómica en poblaciones humanas (Proyecto 1000 Genomas).
   - CEU (Utah residents (CEPH) with Northern and Western European ancestry)
   - FIN (Finnish): consanguinidad por historia poblacional, pequeño número de fundadores y muy poca inmigración.
   - PJL (Punjabi, Pakistan)
   - Cálculo de varios estimadores de consanguinidad con el programa *PLINK*
3. Estimación de la depresión consanguínea en una población ejemplo (*pop_ID*).


## Scripts y datos necesarios

Tras conectarte al CESGA, muévete al directorio de trabajo creado para esta sesión práctica.
```
cd $LUSTRE/MGG_GH5/INBREEDING
pwd
```

Los scripts y datos con los que trabajaremos están almacenados en este mismo repositorio. Para clonar el repositorio desde la terminal:
```
git clone https://github.com/noeliaperezp/MGG_GH5_INBREEDING.git .
ls
```


## 1. Estimación de la consanguinidad genómica en poblaciones simuladas

### *1.1. Simulación de poblaciones*

Utilizaremos el script *script_run_SLiM.sh* para simular los tres escenarios en paralelo. Para enviar el trabajo al sistema de colas se utiliza el comando *sbatch*.

```
sbatch script_run_SLiM.sh slim3INPUT_random_mating
sbatch script_run_SLiM.sh slim3INPUT_bottleneck
sbatch script_run_SLiM.sh slim3INPUT_preference_for_kin
```

Comprueba que se han enviado los trabajos y su estado de ejecución.

```
squeue
```

Una vez finalizado (*squeue* no devuelve ningún trabajo), encontrarás en el directorio *data* dos nuevos archivos para cada uno de los escenarios poblacionales:
```
ls data/
```

- *\*.txt*, con información sobre el tamaño poblacional por generación y heterocigosidad.
- *\*.vcf*, archivo VCF (Variant Call Format) con las variantes simuladas.

### *1.2. Filtrado de datos*

Por lo general, los cálculos realizados por *PLINK* para estimar la consanguinidad no tienen en cuenta el desequilibrio de ligamiento (LD). Además, algunos de los estimadores basados en SNP por SNP son sensibles a las variantes raras. Por este motivo, en ocasiones puede ser recomendable realizar algún tipo de filtrado.

Empezaremos el filtrado con la población simulada bajo panmixia (*slim3INPUT_random_mating*).

Carga los módulos necesarios:
```
module load cesga/2020 gcc/system plink/2.00a2.3
```

Crea un nuevo directorio y muévete a él:
```
mkdir data_pruned
cd data_pruned/
```

Retén variantes bialélicas (SNPs) y genera los archivos con formato plink (.bed,.fam,.bim):
```
plink2 --vcf ../data/slim3INPUT_random_mating.vcf --make-bed --max-alleles 2 --snps-only --out slim3INPUT_random_mating.SNPs
```

Asigna nuevos identificadores a variantes sin nombre:
```
plink2 --bfile slim3INPUT_random_mating.SNPs --set-missing-var-ids @:# --make-bed --out slim3INPUT_random_mating.SNPid
rm slim3INPUT_random_mating.SNPs*
```

Elimina SNPs duplicados:
```
plink2 --bfile slim3INPUT_random_mating.SNPid --rm-dup force-first --make-bed --out slim3INPUT_random_mating.dupfilt
rm slim3INPUT_random_mating.SNPid*
```

Filtra variantes y muestras con >10% missing calls y maf < 0.001:
```
plink2 --bfile slim3INPUT_random_mating.dupfilt --geno 0.1 --mind 0.1 --maf 0.001 --make-bed --out slim3INPUT_random_mating.filtered
rm slim3INPUT_random_mating.dupfilt*
```

Filtra variantes con fuerte desequilibrio de ligamiento (LD):
```
plink2 --bfile slim3INPUT_random_mating.filtered --indep-pairwise 50 5 0.9 --out slim3INPUT_random_mating.plink2
plink2 --bfile slim3INPUT_random_mating.filtered --extract slim3INPUT_random_mating.plink2.prune.in --make-bed --out slim3INPUT_random_mating_pruned
rm slim3INPUT_random_mating.filtered*
rm slim3INPUT_random_mating.plink2*
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/INBREEDING”):
```
cd ..
pwd
```

Para filtrar las poblaciones restantes utilizaremos el mismo código desde un script. Antes de lanzar los siguientes trabajos asegúrate de estar en el directorio correcto, “$LUSTRE/MGG_GH5/INBREEDING”.

```
sbatch script_data_filtering.sh slim3INPUT_bottleneck
sbatch script_data_filtering.sh slim3INPUT_preference_for_kin
```

Verás que todos los datos filtrados se encuentran en el directorio *data_pruned*.
```
ls data_pruned/
```

### *1.3. Estimación de la consanguinidad genómica*

En *PLINK* existen tres comandos que permiten obtener diversos estimadores de consanguinidad. 

El comando *--het* calcula los recuentos de genotipos homocigotos observados y esperados para cada individuo y estima la consanguinidad basada en la desviación de las proporciones de Hardy-Weinberg (*F<sub>LH</sub>*; Li y Horvitz, 1953). El comando *--ibc* calcula diversos estimadores basados en la correlación entre gametos, incluyendo el coeficiente *F* de Yang. 

El comando *--homozyg* infiere consanguinidad a partir de los ROH. Un aspecto destacable de los ROH es que su longitud esperada sigue una distribución exponencial con media 1/(2G) Morgans, donde G es el número de generaciones hasta el ancestro común del que se derivó el segmento. Este aspecto puede utilizarse para distinguir parentesco reciente y lejano. Por ello, utilizaremos tres límites para el tamaño mínimo de los ROH (0.1Mb y 500 generaciones ancestrales, 1Mb y 50 generaciones ancestrales, y 5Mb y 10 generaciones ancestrales, asumiendo una relación 1Mb ≈ 1cM).

Empezaremos con la población simulada bajo panmixia (*slim3INPUT_random_mating_pruned*).

Carga los módulos necesarios:
```
module load cesga/2020 plink/1.9b5
```

Crea un nuevo directorio y muévete a él:
```
mkdir -p results_inbreeding/slim3INPUT_random_mating_pruned
cd results_inbreeding/slim3INPUT_random_mating_pruned/
```

#### --het (F<sub>LH</sub>, F<sub>hom</sub>)
```
plink --bfile ../../data_pruned/slim3INPUT_random_mating_pruned --het --out slim3INPUT_random_mating_pruned.Fhet

sed '1d' slim3INPUT_random_mating_pruned.Fhet.het > temp1       # remove headers
awk '{print ($3/$5)}' temp1 > temp2                             # calculate proportion of homozygous SNPs
sed -i '1i Fhom' temp2                                          # add header
paste slim3INPUT_random_mating_pruned.Fhet.het temp2 > temp3    # add new column 'Fhom'
rm temp1 temp2 slim3INPUT_random_mating_pruned.Fhet.het         # delete unnecessary files
mv temp3 slim3INPUT_random_mating_pruned.Fhet.het               # rename the final output file
```

Encabezados del output **.het*
- FID: ID de la familia
- IID: ID del individuo
- O(HOM): Número observado de homocigotos
- E(HOM): Número esperado de homocigotos
- N(NM): Número de genotipos
- F: Estima del coeficiente de consanguinidad F
- (Añadido) Fhom: Proporción de SNP homocigotos

#### --ibc (F<sub>hatI</sub>, F<sub>hatII</sub>, F<sub>hatIII</sub>)
```
plink --bfile ../../data_pruned/slim3INPUT_random_mating_pruned --ibc --out slim3INPUT_random_mating_pruned.Fibc

column -t slim3INPUT_random_mating_pruned.Fibc.ibc > tmp
rm slim3INPUT_random_mating_pruned.Fibc.ibc
mv tmp slim3INPUT_random_mating_pruned.Fibc.ibc
```
Encabezados del output **.ibc*
- FID: ID de la familia
- IID: ID del individuo
- NOMISS: Número de genotipos
- Fhat1, Fhat2, Fhat3 (*F<sub>Yan</sub>*): Estimas del coeficiente de consanguinidad *F*

#### --homozyg (F<sub>ROH</sub>, ROH > 0.1 Mb)
```
plink --bfile ../../data_pruned/slim3INPUT_random_mating_pruned --homozyg-kb 100 --homozyg --out slim3INPUT_random_mating_pruned.FROH01

sed '1d' slim3INPUT_random_mating_pruned.FROH01.hom.indiv > temp1     # remove headers
awk '{print ($5/100000)}' temp1 > temp2                               # calculate FROH (genome length 100Mb = 100000Kb)
sed -i '1i FROH01' temp2                                              # add header
paste slim3INPUT_random_mating_pruned.FROH01.hom.indiv temp2 > temp3  # add new column 'FROH'
rm temp1 temp2 slim3INPUT_random_mating_pruned.FROH01.hom.indiv       # delete unnecessary files
mv temp3 slim3INPUT_random_mating_pruned.FROH01.hom.indiv             # rename the final output file
```

Encabezados del output **.hom*
- FID: ID de la familia
- IID: ID del individuo
- PHE: fenotipo
- CHR: Cromosoma
- SNP1: SNP al inicio del segmento
- SNP2: SNP al final del segmento
- POS1: Posición física (bp) del SNP1
- POS2: Posición física (bp) del SNP2
- KB: Longitud del segmento (kb)
- NSNP: Número de SNPs en ROH
- DENSITY: Densidad media de SNPs (1 SNP per kb)
- PHOM: Proporción de sitios homocigotos
- PHET: Proporción de sitios heterocigotos

Encabezados del output **.hom.indiv*
- FID: ID de la familia
- IID: ID del individuo
- PHE: fenotipo
- NSEG: Número de ROHs por individuo
- KB: longitud total en ROH
- KBAVG: Tamaño medio de los ROH
- (añadido) FROH: estima del coeficiente de consanguinidad

#### --homozyg (F<sub>ROH</sub>, ROH > 1 Mb)
```
plink --bfile ../../data_pruned/slim3INPUT_random_mating_pruned --homozyg-kb 1000 --homozyg --out slim3INPUT_random_mating_pruned.FROH1

sed '1d' slim3INPUT_random_mating_pruned.FROH1.hom.indiv > temp1     # remove headers
awk '{print ($5/100000)}' temp1 > temp2                              # calculate FROH (genome length 100Mb = 100000Kb)
sed -i '1i FROH1' temp2                                              # add header
paste slim3INPUT_random_mating_pruned.FROH1.hom.indiv temp2 > temp3  # add new column 'FROH'
rm temp1 temp2 slim3INPUT_random_mating_pruned.FROH1.hom.indiv       # delete unnecessary files
mv temp3 slim3INPUT_random_mating_pruned.FROH1.hom.indiv             # rename the final output file
```

#### --homozyg (F<sub>ROH</sub>, ROH > 5 Mb)
```
plink --bfile ../../data_pruned/slim3INPUT_random_mating_pruned --homozyg-kb 5000 --homozyg --out slim3INPUT_random_mating_pruned.FROH5

sed '1d' slim3INPUT_random_mating_pruned.FROH5.hom.indiv > temp1     # remove headers
awk '{print ($5/100000)}' temp1 > temp2                              # calculate FROH (genome length 100Mb = 100000Kb)
sed -i '1i FROH5' temp2                                              # add header
paste slim3INPUT_random_mating_pruned.FROH5.hom.indiv temp2 > temp3  # add new column 'FROH'
rm temp1 temp2 slim3INPUT_random_mating_pruned.FROH5.hom.indiv       # delete unnecessary files
mv temp3 slim3INPUT_random_mating_pruned.FROH5.hom.indiv             # rename the final output file
```

Ahora que tenemos todas las estimas inferidas, generaremos una tabla resumen con todos los valores:
```
awk '{print $2 "\t" $7}' slim3INPUT_random_mating_pruned.Fhet.het > temp1     # FHOM
awk '{print $6}' slim3INPUT_random_mating_pruned.Fibc.ibc > temp2             # Fhat3 (or FYan)
awk '{print $7}' slim3INPUT_random_mating_pruned.FROH01.hom.indiv > temp3     # FROH01
awk '{print $7}' slim3INPUT_random_mating_pruned.FROH1.hom.indiv > temp4      # FROH1
awk '{print $7}' slim3INPUT_random_mating_pruned.FROH5.hom.indiv > temp5      # FROH5
paste temp1 temp2 temp3 temp4 temp5 | column -t > slim3INPUT_random_mating_pruned.inbreeding
rm temp1 temp2 temp3 temp4 temp5
```

Retorna al directorio principal (“$LUSTRE/MGG_GH5/INBREEDING”):
```
cd ../..
pwd
```

Para estimar la consanguinidad en las demás poblaciones simuladas utilizaremos los mismos comandos desde un script.
```
sbatch script_estimate_inbreeding.sh slim3INPUT_bottleneck_pruned
sbatch script_estimate_inbreeding.sh slim3INPUT_preference_for_kin_pruned
```

Consanguinidad promedio en las tres poblaciones simuladas:
```
# Panmixia
awk 'NR==1 {print; next} {for(i=1;i<=NF;i++) sum[i]+=$i} END {for(i=1;i<=NF;i++) printf "%s%.4f", (i==1?"":" "), sum[i]/(NR-1); print ""}' results_inbreeding/slim3INPUT_random_mating_pruned/slim3INPUT_random_mating_pruned.inbreeding | column -t

# Cuello de botella
awk 'NR==1 {print; next} {for(i=1;i<=NF;i++) sum[i]+=$i} END {for(i=1;i<=NF;i++) printf "%s%.4f", (i==1?"":" "), sum[i]/(NR-1); print ""}' results_inbreeding/slim3INPUT_bottleneck_pruned/slim3INPUT_bottleneck_pruned.inbreeding | column -t

# Preferencia por apareamientos entre parientes
awk 'NR==1 {print; next} {for(i=1;i<=NF;i++) sum[i]+=$i} END {for(i=1;i<=NF;i++) printf "%s%.4f", (i==1?"":" "), sum[i]/(NR-1); print ""}' results_inbreeding/slim3INPUT_preference_for_kin_pruned/slim3INPUT_preference_for_kin_pruned.inbreeding | column -t
```

## 2. Estimación de la consanguinidad genómica en poblaciones humanas (Proyecto 1000 Genomas)

Estima la consanguinidad para cada una de las poblaciones humanas, CEU, FIN, y PJL, siguiendo el mismo procedimiento que con las poblaciones simuladas.

Filtrado de datos: 
```
sbatch script_data_filtering.sh chr22_CEU
sbatch script_data_filtering.sh chr22_FIN
sbatch script_data_filtering.sh chr22_PJL
```

Cálculo de la consanguinidad (teniendo en cuenta que la longitud del cromosoma 22 es de aproximadamente 51Mb):
```
sbatch script_estimate_inbreeding.sh chr22_CEU_pruned
sbatch script_estimate_inbreeding.sh chr22_FIN_pruned
sbatch script_estimate_inbreeding.sh chr22_PJL_pruned
```

Consanguinidad promedio por población:
```
# CEU
awk 'NR==1 {print; next} {for(i=1;i<=NF;i++) sum[i]+=$i} END {for(i=1;i<=NF;i++) printf "%s%.4f", (i==1?"":" "), sum[i]/(NR-1); print ""}' results_inbreeding/chr22_CEU_pruned/chr22_CEU_pruned.inbreeding | column -t

# FIN
awk 'NR==1 {print; next} {for(i=1;i<=NF;i++) sum[i]+=$i} END {for(i=1;i<=NF;i++) printf "%s%.4f", (i==1?"":" "), sum[i]/(NR-1); print ""}' results_inbreeding/chr22_FIN_pruned/chr22_FIN_pruned.inbreeding | column -t

# PJL
awk 'NR==1 {print; next} {for(i=1;i<=NF;i++) sum[i]+=$i} END {for(i=1;i<=NF;i++) printf "%s%.4f", (i==1?"":" "), sum[i]/(NR-1); print ""}' results_inbreeding/chr22_PJL_pruned/chr22_PJL_pruned.inbreeding | column -t
```


## 3. Estimación de la depresión consanguínea en una población ejemplo (*pop_ID*).

Para poder estimar la depresión consanguínea necesitamos, además de los coeficientes de consanguinidad, fenotipos para algún caracter de la eficacia. Estos datos se proporcionan en el ejemplo *pop_ID*.

Filtrado de datos
```
sbatch script_data_filtering.sh pop_ID
```

Cálculo de la consanguinidad
```
sbatch script_estimate_inbreeding.sh pop_ID_pruned
```

La pendiente de la regresión lineal proporciona una estima del cambio en la media por unidad de *F*. 


### *Transferencia de ficheros del CESGA a tu sistema local*

Representaremos los resultados en RStudio. Para ello necesitarás transferir los resultados del cálculo de la consanguinidad molecular a tu sistema local. También necesitarás el archivo **.phe* con los fenotipos de cada individuo. Para transferir estos archivos puedes usar el programa WinSCP, o alternativamente usando el comando *scp* desde tu máquina local. Recuerda cambiar:

- *username* por tu nombre de usuario
- */path* por la ruta que te devuelve *pwd* en el CESGA

```
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/results_inbreeding/pop_ID_pruned/pop_ID_pruned.inbreeding .
scp -o MACs=hmac-sha2-512 username@ft3.cesga.es:/path/data/pop_ID.phe .
```

Una vez transferidos los resultados, abre RStudio y asegúrate de que te encuentras en la carpeta que contiene los resultados. A continuación copia y ejecuta las líneas de código que encontrarás en el script inbreeding_figures.R.

