# HRV-WA
HRV-WA algorithm, based on Frequency Möbius  Modulated (FMM) models [1], analyses 24-h HRV rhythms from RR data. HRV-WA predicts HRV and identifies Nocturnal and Diurnal waves involved in HRV analysis, see [citation] for details.

# How to use 
HRV-WA code is implemented in R language programming and is easy to use. 
The code provided in this repository replicates the steps described in [citation] where HRV-WA was described.

Open the R script named ``runHRV-WA.R`` including an example to conduct the methodology. 

First, source R functions saved in ``codeHRVWA.R`` must be loaded. 

````
source("codeHRVWA.R")
````

Then, load the file ``dataExample.RData``. It contains an example of 24-h RR data across 5-min time intervals that serves as input.

INPUTS: 
````
load(file="dataExample.RData")
````
HRV-WA algorithm is applied on the data using the function ``HRV_WA``.

RUN:
````
HRV_WA(data=dataExample)
````

A matrix with the FMM parameter estimations for Principal and Secondary FMM wave(s) is given. 

OUTPUT:
````
     Comp         M         A    Alpha     Beta      Omega       t_U      t_L       R2_m
Nocturnal    1 0.7653564 0.2102229 4.335471 4.168209 0.57085754 0.5712119 2.777274 0.79368260
Diurnal    3 0.7653564 0.1114734 2.248823 4.381641 0.05958329 5.3053883 5.556942 0.04729768
````

In addition, plots with the HRV prediction and wave decomposition are also provided. 
![Pac_2_intro](https://user-images.githubusercontent.com/24298539/214274483-c23af48c-2ca8-48aa-a4bb-677f4abedbf1.jpg)

# References
[1] Rueda, C., Larriba, Y. & Peddada, S.D. Frequency Modulated Möbius Model Accurately Predicts Rhythmic Signals in Biological and Physical Sciences. Sci Rep 9, 18701 (2019). https://doi.org/10.1038/s41598-019-54569-1

[citation]
