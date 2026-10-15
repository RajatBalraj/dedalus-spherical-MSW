How this module works:

1. 
we start at dedalus_grid_setup.py that is called in every subsequent module
this is where the units and the grid, basis, coords, variables get defined
so unit constants remain consistent across modules
Gives a class object whose attributes are assigned to each of the required structures in each module

2.
Then comes the IC file. For now it is gaussian_height_qsat_IC, but MAY BE CHANGED to something else also. This will create the initial height/Q perturbation and the saturation moisture field. I suppose later we'll have to make two separate files for the perturb and the q_sat

This is necessary because these variables do not get created in the file no. 1, nor in the actual run file. They do not get saved in the spinup checkpoints either. Hence it will be called in both the subsequent modules

the random_height_qsat_IC is different. It saves the h and q_sat into a startup file that has to be loaded in using solver.state. This is because the random forcing is generated without using pure dedalus operations. So computing it in parallel ends up in shape mismatches and what not. There should be a better waybut this is the best fix i've figured out for now.

3. spinup.py
Executes the spinup and saves a dry and moist part. Calls IC file. If you change the IC file you need to modify this.
later perhaps we can shift that decision to the terminal but we'll see

Its inputs take IC_type of random or gaussian. In the random case the dlon and lon_center cases are useless and any value can be entered in there

4. actual_run.py
Also configured to the gaussian IC but agnostic to the details. Has two options, spinup 1 or 0, this will determine if the full equations run is done with the spun up case or with only the perturb IC with everything else set to zero respectively




BATCH file

taskset -c 0-23 mpirun -n 24 --bind-to none python full_run_batch_script.py  --data_folder=random_forcing_dlat30 --run_days=300 --dlat=30 --lat_center=0 --spinup=0 --IC_type=random --lon_center=0 --dlon=0 > random_perturb.out


