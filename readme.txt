To verify the mass conservation:
A simple H permeation benchmark is performed, consisting of 2 phases:
	1. Charging phase:
		- During: 200s
		- Constant temperature and flux: 293 K
		- Constant flux: 8e11 H/mm2.s
		
	2. Heating phase: to observe the desorption behaviour and verify mass conservation
		- During: 800s
		- Linear temperature ramp: from 293K to 793K with heating rate of 0.625K/s
		- Zero H flux on all boundaries.

To reproduce the result:
	1. Run input file " Job-1.inp " with the user subroutine file "UMATHT_transient_void.f ".
	2. Post-processing for the evolution of H retention: run python script " postProc.py " in Abaqus.
	3. Plot the results: run python script " plot.py " in python environment.
