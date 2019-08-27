# Computationally inexpensive identification of noninformative model parameters by sequential screening
*by Matthias Cuntz (INRA Nancy, France) and Juliane Mai (University of Waterloo, Canada) et al.*

## Abstract
Environmental models tend to require increasing computational time and resources as physical process descriptions are improved or new descriptions are incorporated. Many-query applications such as sensitivity analysis or model calibration usually require a large number of model evaluations leading to high computational demand. This often limits the feasibility of rigorous analyses. Here we present a fully automated sequential screening method that selects only informative parameters for a given model output. The method is called **Efficient Elementary Effects** and requires a number of model evaluations that is approximately 10 times the number of model parameters. It was tested using the mesoscale hydrologic model mHM in three hydrologically unique European river catchments. It identified around 20 informative parameters out of 52, with different informative parameters in each catchment. The screening method was evaluated with subsequent analyses using all 52 as well as only the informative parameters. Subsequent Sobol’s global sensitivity analysis led to almost identical results yet required 40% fewer model evaluations after screening. mHM was calibrated with all and with only informative parameters in the three catchments. Model performances for daily discharge were equally high in both cases with Nash-Sutcliffe efficiencies above 0.82. Calibration using only the informative parameters needed just one third of the number of model evaluations. The universality of the sequential screening method was demonstrated using several general test functions from the literature. We therefore recommend the use of the computationally inexpensive sequential screening method prior to rigorous analyses on complex environmental models. The full paper can be found [here](https://doi.org/10.1002/2015WR016907).

## Step-by-Step Tutorial
The step-by-step tutorial describes all the steps to apply the method of **Efficient Elementary Effects (EEE)** in order to identify the set of informative parameters of a model in a fully automated manner. Details can be found [here](https://github.com/julemai/EEE/wiki/Step-by-Step-Tutorial).

## Examples
To be filled.

## Citation
Cuntz, M., & Mai, J. et al. (2015).<br>
Computationally inexpensive identification of noninformative model parameters by sequential screening.<br>
*Water Resources Research*, 51, 6417–6441.<br>
https://doi.org/10.1002/2015WR016907.

