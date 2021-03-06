[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3576613.svg)](https://doi.org/10.5281/zenodo.3576613)

# Heterogeneous_SOMDynamics

This repository contains Matlab R2018b scripts to simulate soil organic matter (SOM) dynamics in a heterogeneous soil system and is a part of our publication [Chakrawal et al. 2020](https://www.geosci-model-dev.net/13/1399/2020/gmd-13-1399-2020.html). **main_XXX.m** files are in Scenario1_SteadyStateIC and Scenario2_transientIC folders and solve the mass balance equations of SOM compartments. In the folder Scenario1_SteadyStateIC, a soil system initially on average in steady state conditions is considered; in the folder Scenario2_transientIC, the responses to initial conditions perturbed with respect to far from the steady state are studied. In both scenarios, Michaelis-Menten and Multiplicative kinetics for decomposition are compared.
![Fig1](https://github.com/ArjunChakrawal/Heterogeneous_SOMDynamics/blob/master/GMD_fig1.png)

The folder structure is as follows:

- SOM_heteogeneous
    - Scenario1_SteadyStateIC
        - results
        - Spatial_field
    - Scenario2_transientIC
        - MM_Inv_transient
		- MM_transient
		- Mult_transport
        - Spatial_field
    - Third_party_scripts


To run the code, download the zip file of the repository and unzip it to your favourite location, and afterwards follow the instructions below. Here we show the example of Scenario1_SteadyStateIC.
<br/>
### Step1: Spatial field generation
Run *Scenario1_SteadyStateIC/Spatial_field/SS_field_generator.m* to generate heterogeneous fields of substrate and microbial C. While running this code, a user interface opens to select the **Heterogeneous_SOMDynamics** folder. IMPORTANT: For reproducibility of results, the provided ksmmh1.mat, kmh1.mat and ksm1.mat files, which contain spatially heterogeneous maps for kinetics parameters, should be used. By default save commands to save the heterogeneous fields files are commented. To create new heterogeneous fields uncomment the save commands in SS_field_generator.m file.
### Step2: SOM dynamics in heterogeneous domain
Run *Scenario1_SteadyStateIC/main_Mult.m* to simulate SOM dynamics in the heterogeneous domain created in Step 1. While running this code, a user interface opens to select the **Heterogeneous_SOMDynamics** folder.
### Step3: Post-processing
To plot the results scripts are provided coresspondig folders named as **Fig_X.m**. While running this code, a user interface opens to select the **Heterogeneous_SOMDynamics** folder.

<br/>
Some third party scripts (from MATLAB Central File Exchange) are used in delevoping this model and already included in this repository. These scripts are redistributable under their respective licenses which can be found in the following links.

* [parfor_progressbar](https://www.mathworks.com/matlabcentral/fileexchange/53773-parfor_progressbar), by Daniel Terry 2016. 
* [align_Ylabels](https://se.mathworks.com/matlabcentral/fileexchange/41701-y-labels-alignment-in-subplots?focused=3788739&tab=function), by Giuliano Bernardi 2013. 
* [spatialPattern](https://se.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data), by Jon Yearsley 2016.
* [linspecer]( https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap), by  Jonathan C. Lansey 2015.
* [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig), by Yair Altman 2018
* [tight_subplot](https://se.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w), by  Pekka Kumpulainen 2016.
* [ds2nfu](https://se.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion), by Michelle Hirsch 2015.

## Citation
Please use the following citation [GMD](https://www.geosci-model-dev.net/13/1399/2020/gmd-13-1399-2020.html).
```
@Article{gmd-13-1399-2020,
AUTHOR = {Chakrawal, A. and Herrmann, A. M. and Koestel, J. and Jarsj\"o, J. and Nunan, N. and K\"atterer, T. and Manzoni, S.},
TITLE = {Dynamic upscaling of decomposition kinetics for carbon cycling models},
JOURNAL = {Geoscientific Model Development},
VOLUME = {13},
YEAR = {2020},
NUMBER = {3},
PAGES = {1399--1429},
URL = {https://www.geosci-model-dev.net/13/1399/2020/},
DOI = {10.5194/gmd-13-1399-2020}
}
}
```
