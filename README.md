# Heterogeneous_SOMDynamics

This repository contains Matlab R2018b scripts to simulate soil organic matter (SOM) dynamics in a heterogeneous soil system. **main_XXX.m** files are in Scenario1_SteadyStateIC and Scenario2_transientIC folders and solve the mass balance equations of SOM compartments. In the folder Scenario1_SteadyStateIC, a soil system initially on average in steady state conditions is considered; in the folder Scenario2_transientIC, the responses to initial conditions perturbed with respect to far from the steady state are studied. In both scenarios, Michaelis-Menten and Multiplicative kinetics for decomposition are compared.

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
* [spatialPattern](https://se.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data), by Jon Yearsley 2016.
* [linspecer]( https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap), by  Jonathan C. Lansey 2015.
* [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig), by Yair Altman 2018
* [tight_subplot](https://se.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w), by  Pekka Kumpulainen 2016.
* [ds2nfu](https://se.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion), by Michelle Hirsch 2015.
