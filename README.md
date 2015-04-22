# rocket-ode
Simulation of water bottle rocket physics as ODEs

I created this model for the final portion of ME495 W15.

The following files are contained here:
* odel_sim.m			*calls main simulation function modelfun. Requires eventLA.m and eventZ0.m*
* modelfun.m			*main function for simulation*
* launch_val.m			*calls launch validation function valfun. Requires eventLA.m and eventZ0.m*
* valfun.m				*main function for launch validation*
* simplefun.m			*simplified version of modelfun.m*
* eventLA.m	event 		*event function for rocket leaving the sting*
* eventZ0.m	event 		*event function for rocket hitting the ground*

Instructions to run a simulation:
1. Put model_sim.m, modelfun.m, eventLA.m, and eventZ0.m in the same folder (replace modelfun.m with simplefun.m if desired)
2. Open model_sim.m and enter desired initial conditions
3. (Some constant conditions, such as wind speed, are in modelfun.m. 4. Un-comment desired plot outputs
5. Run model_sim.m
6. MATLAB will output final distance and any uncommented plots

Instructions to run a launch validation:
1. Put launch_val.m, valfun.m, eventLA.m, and eventZ0.m in the same folder
2. Run launch_val.m
3. MATLAB will plot the velocity vs time and position vs time of the model and the sample data

