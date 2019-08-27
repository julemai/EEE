# Computationally inexpensive identification of noninformative model parameters by sequential screening
*by Matthias Cuntz (INRA Nancy, France) and Juliane Mai (University of Waterloo, Canada) et al.*

## Abstract
Environmental models tend to require increasing computational time and resources as physical process descriptions are improved or new descriptions are incorporated. Many-query applications such as sensitivity analysis or model calibration usually require a large number of model evaluations leading to high computational demand. This often limits the feasibility of rigorous analyses. Here we present a fully automated sequential screening method that selects only informative parameters for a given model output. The method is called **Efficient Elementary Effects (EEE)** and requires a number of model evaluations that is approximately 10 times the number of model parameters. It was tested using the mesoscale hydrologic model mHM in three hydrologically unique European river catchments. It identified around 20 informative parameters out of 52. The universality of the sequential screening method was demonstrated using several general test functions from the literature. The full paper can be found [here](https://doi.org/10.1002/2015WR016907).

## Step-by-Step Tutorial
The step-by-step tutorial describes all the steps to apply the method of **Efficient Elementary Effects (EEE)** in order to identify the set of informative parameters of a model in a fully automated manner. Details can be found [here](https://github.com/julemai/EEE/wiki/Step-by-Step-Tutorial).

## Examples
We provide a few example workflows on how to use the provided codes in order to obtain the non informative parameters of some models using EEE. Details can be found [here](https://github.com/julemai/EEE/wiki/Examples).

## Citation
M Cuntz & J Mai et al. (2015).<br>
Computationally inexpensive identification of noninformative model parameters by sequential screening.<br>
*Water Resources Research*, 51, 6417–6441.<br>
https://doi.org/10.1002/2015WR016907.
