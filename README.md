# HRV-WA
HRV-WA algorithm, based on Frequency Möbius  Modulated (FMM) models, analyses 24-h HRV rhythms directly from RR data. HRV-WA not only predicts HRV, but also identifies and assigns Direct, and at times, Guided wave(s) involved in the analysis.

# How to use 
HRV-WA is achieved in R and is easy to use. 
The code provided in this GitHub replicated the four steps described in [citation].

Run the R script named runHRV-WA.R to ,conduct the methodology.
The file dataExample.RData is loaded on the R script and serves as example of 24-h RR data across 5-min time intervals.

INPUTS: 
  - Vector of RR data (data).

RUN:
````
load(file="dataExample.RData")
HRV_WA(data=dataExample)
````

OUTPUT:
  - Matrix with the FMM parameter estimations for Direct (and Guided) FMM wave(s) assigned.
  ````
     Comp         M         A    Alpha     Beta      Omega       t_U      t_L       R2_m
Direct    1 0.7653564 0.2102229 4.335471 4.168209 0.57085754 0.5712119 2.777274 0.79368260
Guided    3 0.7653564 0.1114734 2.248823 4.381641 0.05958329 5.3053883 5.556942 0.04729768
````
  - Plots with the HRV prediction and wave decomposition. 
![image](https://user-images.githubusercontent.com/24298539/209979899-31c967de-408a-44e4-b1ee-60a0c9425fe4.png)
