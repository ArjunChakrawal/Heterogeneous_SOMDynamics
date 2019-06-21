# Heterogeneous_SOMDynamics


This repository contain matlab script to run soil organic matter (SOM) dynamics in a hetergeneous soil system.
main files are in Scenario1_SteadyStateIC and Scenario2_transientIC folders for simulating SOM dynamic using Michaelis-Menten
and Multiplicative kinetics.

Folder structure is following

- SOM_heteogeneous
    - Scenario1_SteadyStateIC
        - results
        - Spatial_field
    - Scenario2_transientIC
        - results
        - Spatial_field

Code must be run in the following order. Here we use the example of Scenario1_SteadyStateIC.
### Step1: Spatial field generation
Run *\Scenario1_SteadyStateIC\Spatial_field\SS_field_generator.m* to generate heterogeneous fields of substrate and microbial C.
IMPORTANT: For reproducebility of results, provided *ksmmh1.mat*, *kmh1.mat* and *ksm1.mat* files should be used.
### Step2: SOM dynamics in heterogeneous domain
Run *\Scenario1_SteadyStateIC\main_Mult.m* file to simulate SOM dynamics in a hetergeneous soil.
### Step3: Post-processing
Some scripts are provide for post-processing of results in the 
Results are stored in the *\Scenario1_SteadyStateIC\results\* folder. 

<br/>
There are some third party scripts used in delevoping this model and already included in this repository.

* [parfor_progressbar](https://www.mathworks.com/matlabcentral/fileexchange/53773-parfor_progressbar)
* [spatialPattern](https://se.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data)
* [linspecer]( https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap)
