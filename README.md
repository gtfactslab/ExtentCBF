# Extent-Compatible Control Barrier Functions

This is code for the extent-compatible control barrier functions framework. This repository contains two implementations of the proposed framework in the Robotarium simulator:

* A sum-of-squares (SOS) based approach
* A sampling-based approach implemented on the Robotarium

A classical CBF implementation where one has to shrink the safe set is also included for refernece.


## To execute the code, please follow the instructions below:

* Enter the Robotarium_Extent_Experiment folder  
* Run the init.m file
* In implementation/main.m: Modify the variable cont to be either 'extent', 'sos', 'point' to choose which controller that should be used
* To run the simulation, just run implementation/main.m

NOTE: The SOS controller code requires SOSTOOLS with solver SDPT3. It has only been varified to work with SOSTOOLS v 3.03, SDPT3 v 4.0 and MATLAB R2019b.