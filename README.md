# Heterogeneous_SOMDynamics

This repository contains Matlab 2018b scripts to simulate soil organic matter (SOM) dynamics in a heterogeneous soil system. **main_XXX.m** files are in Scenario1_SteadyStateIC and Scenario2_transientIC folders and solve the mass balance equations of SOM compartments. In the folder Scenario1_SteadyStateIC, a soil system initially on average in steady state conditions is considered; in the folder Scenario2_transientIC, the responses to initial conditions perturbed with respect to the steady state are studied. In both scenarios, Michaelis-Menten and Multiplicative kinetics for decomposition are compared.

The folder structure is as follows:

- SOM_heteogeneous
    - Scenario1_SteadyStateIC
        - results
        - Spatial_field
    - Scenario2_transientIC
        - results
        - Spatial_field
    - Third_party_scripts


To run the code, download the zip file of the repository and unzip it to your favourite location, and afterwards follow the instructions below. Here we show the example of Scenario1_SteadyStateIC.
<br/>
### Step1: Spatial field generation
Run *Scenario1_SteadyStateIC\Spatial_field\SS_field_generator.m* to generate heterogeneous fields of substrate and microbial C. While running this code, a user interface opens to select the **Heterogeneous_SOMDynamics** folder. IMPORTANT: For reproducibility of results, the provided ksmmh1.mat, kmh1.mat and ksm1.mat files, which contain spatially heterogeneous maps for kinetics parameters, should be used.
### Step2: SOM dynamics in heterogeneous domain
Run *Scenario1_SteadyStateIC\main_Mult.m* to simulate SOM dynamics in the heterogeneous domain created in Step 1. While running this code, a user interface opens to select the **Heterogeneous_SOMDynamics** folder.
### Step3: Post-processing
Some scripts are provide for the post-processing of results in the *Scenario1_SteadyStateIC\results* folder. While running this code, a user interface opens to select the **Heterogeneous_SOMDynamics** folder.

<br/>
There are some third party scripts used in delevoping this model and already included in this repository. All these scripts are redistributabel but under thier respective licenses provide in the links.

* [parfor_progressbar](https://www.mathworks.com/matlabcentral/fileexchange/53773-parfor_progressbar)
* [spatialPattern](https://se.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data)
* [linspecer]( https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap)
* [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
* [tight_subplot](https://se.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
* [ds2nfu](https://se.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion)
