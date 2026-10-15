1. INITIALISING GRID
we start at dedalus_grid_setup.py that is called in every subsequent module
this is where the units and the grid, basis, coords, variables get defined so unit constants remain consistent across modules
Gives a class object whose attributes are assigned to each of the required structures in each module

2. CREATING INITIAL CONDITION
(a) Then comes the IC file. For now it is gaussian_height_qsat_IC and random_height_qsat_IC to something else also. This will create the initial height perturbation and the saturation moisture field. I suppose later we'll have to make two separate files for the perturb and the q_sat.

(b) Random height forcing is defined by summing up first n(given) spherical wavenumbers (large scale)in random linear combination and adding on a decay function so the disturbance is limited to the given latitudes 

(c) This file is necessary because these variables do not get created in the file no. 1, nor in the actual run file. They do not get saved in the spinup checkpoints either. Hence it will be called in both the subsequent modules.

(d) Due to dedalus having its own parallel computing workflow, the random IC and gaussian IC need to be generated differently. Gaussian can be created with a paralellized process but random IC requires creating and saviing the IC first and then calling it via dedalus handlers later in subsequent processes.

(e) The random_height_qsat_IC saves the h and q_sat into a startup file that has to be loaded in using solver.state. This is because the random forcing is generated without using pure dedalus operations. So computing it in parallel ends up in shape mismatches and other issues. There could be a better way but this is the best fix I've figured out for now.

3. SPINUP
Executes the spinup and saves a dry and moist part. Calls IC file. If you change the IC file you need to modify this.
Its inputs take IC_type of random or gaussian. In the random case the dlon and lon_center cases are useless and any value can be entered in there
Spunup up state is arrived at by holding the initally applied height perturbation constant and solving for the other variables. This gives us a quasi-balanced initial state

5. FINAL RUN (actual_run.py)
Has two options, spinup 1 or 0, this will determine if the full equations run is done with the spun up case or with only the perturb height IC with everything else set to zero respectively.

6. BATCH file

Options to initial with gaussian or random height perturbation. Option to execute spinup or not. Executes above files in order as required.
Example terminal command: taskset -c 0-15 mpirun -n 16 --bind-to none python full_run_batch_script.py  --data_folder=random_forcing_dlat30 --run_days=300 --dlat=30 --lat_center=0 --spinup=0 --IC_type=random --lon_center=0 --dlon=0 > random_perturb.out


