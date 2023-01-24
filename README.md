# HRV-WA
HRV-WA algorithm, based on Frequency Möbius  Modulated (FMM) models [1], analyses 24-h HRV rhythms from RR data. HRV-WA predicts HRV and identifies Principal and Secondary waves involved in HRV analysis, see [citation] for details.

# How to use 
HRV-WA is achieved in R and is easy to use. 
The code provided in this repositury replicates the steps described in [citation] where HRV-WA was described.

Run the R script named runHRV-WA.R to conduct the methodology.
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

# References
[1] Rueda, C., Larriba, Y. & Peddada, S.D. Frequency Modulated Möbius Model Accurately Predicts Rhythmic Signals in Biological and Physical Sciences. Sci Rep 9, 18701 (2019). https://doi.org/10.1038/s41598-019-54569-1
[citation]
