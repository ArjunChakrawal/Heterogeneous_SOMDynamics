# Heterogeneous_SOMDynamics


This repository contains matlab scripts to simulate soil organic matter (SOM) dynamics in a hetergeneous soil system.
**main_XXX.m** files are in Scenario1_SteadyStateIC and Scenario2_transientIC folders for solving mass balance equations using Michaelis-Menten and Multiplicative kinetics.

Folder structure is following

- SOM_heteogeneous
    - Scenario1_SteadyStateIC
        - results
        - Spatial_field
    - Scenario2_transientIC
        - results
        - Spatial_field
    - Third_party_scripts


NOTE: Code must be run in the following order. Here we show the example of Scenario1_SteadyStateIC.
<br/>
Download the zip file of the repository and unzip it to favourite location, afterwards follow the instructions below.
### Step1: Spatial field generation
Run *Scenario1_SteadyStateIC\Spatial_field\SS_field_generator.m* to generate heterogeneous fields of substrate and microbial C. While running this code, it opens a UI to select the **Heterogeneous_SOMDynamics** folder.
IMPORTANT: For reproducebility of results, provided *ksmmh1.mat*, *kmh1.mat* and *ksm1.mat* files should be used.
### Step2: SOM dynamics in heterogeneous domain
Run *Scenario1_SteadyStateIC\main_Mult.m* file to simulate SOM dynamics in a hetergeneous soil. While running this code, it opens a UI to select the **Heterogeneous_SOMDynamics** folder.
### Step3: Post-processing
Some scripts are provide for the post-processing of results in the *Scenario1_SteadyStateIC\results\* folder. While running this code, it opens a UI to select the **Heterogeneous_SOMDynamics** folder.

<br/>
There are some third party scripts used in delevoping this model and already included in this repository.

* [parfor_progressbar](https://www.mathworks.com/matlabcentral/fileexchange/53773-parfor_progressbar)
* [spatialPattern](https://se.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data)
* [linspecer]( https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap)
* [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
* [tight_subplot](https://se.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
* [ds2nfu](https://se.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion)
