
# MGG GH5 - ANÁLISIS DE CONSANGUINIDAD
* breve introducción teórica de la consanguinidad y los estimadores. Depresión consanguínea

* breve explicación de la utilidad de las simulaciones

* breve explicación de slim

* breve mención a los datos humanos que se utilizarán. Consanguinidad en poblaciones humanas

* Tabla resumen de F ejemplo resultado de apareamientos entre familiares.


## Scripts necesarios

Tras conectarte al CESGA, muévete al directorio de trabajo creado para esta sesión práctica.
```
cd $LUSTRE/MGG_GH5/INBREEDING
```

Los scripts con los que trabajaremos están almacenados en este mismo repositorio. Para clonar el repositorio desde la terminal:
```
git clone https://github.com/noeliaperezp/MGG_GH5_INBREEDING.git .
ls
```


## Simulación de poblaciones
* breve descripción del script
* descripción de los escenarios (figura de N para mostrar el cuello de botella)

```
sbatch script_run_SLiM.sh slim3INPUT_random_mating
```

```
sbatch script_run_SLiM.sh slim3INPUT_bottleneck
```

```
sbatch script_run_SLiM.sh slim3INPUT_partial_selfing
```

```
sbatch script_run_SLiM.sh slim3INPUT_preference_for_kin
```

```
squeue
```

breve mención de los resultados que salen y formato


## Cálculo de consanguinidad genómica


## Estimación de la depresión consanguínea


## Consanguinidad en poblaciones humanas


