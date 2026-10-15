import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

import sys
# 1. Provide the path to the folder containing your .py file
sys.path.append('/home/rajat/Desktop/Project/codes/random_forcing')

import sphericalForcing

import dedalus_grid_setup
import xarray as xr
dtype = np.float64
import argparse

data_path = '/home/rajat/Desktop/Project/data/'

parser = argparse.ArgumentParser(description="target save folder and spinup IC choice")
parser.add_argument("--data_folder", type=str,required=True,
                    help="folder to which to save all data. Will also take the initial condition from the checkpoint files in this folder")


parser.add_argument("--dlat", type=float, required=True, help="lat extent of gaussian in degrees")
parser.add_argument("--lat_center", type=float, default=0, help="lat center of gaussian in degrees")

args = parser.parse_args()
#################################################################
# Simulation units

ug = dedalus_grid_setup.units_and_grid(300, 0.05, 5000, 0.1) #units and grids #chi and tau fields are just placeholders
#, not relevant to this module
meter = ug.meter
dist = ug.dist
Q = ug.Q #scale of saturation vapor pressure in meters of water in the column (50mm water about how much is found in the equatorial atm)
h = ug.h
q_sat = ug.q_sat
lat_field = ug.lat_field

##########################################################################

# Initial conditions: gaussian bump


if dist.comm.rank==0:
    perturb = sphericalForcing.sphForcing(512, 256, 8, 6371e3, 0, 12, 2, 0)
    '''
    nlons, nlats, ntrunc, rsphere, lmin, lmax, magnitude, correlation

    #ntrunc       # typically 2/3rd of max resolved wavenumber
        magnitude given as 2 meters, #as a rule of thumb, the final forcing amplitude seems to end up 3-5 times the number here. so check before applying to simulation
        correlation is zero - white noise
        
    '''

    perturb_actual = perturb.forcingFn(0)
    my_sphere = perturb.trans

    #forcing_real = my_sphere.spectogrd(perturb_actual) 
    lat_center = np.deg2rad(args.lat_center)
    dlat = np.deg2rad(args.dlat)

    lats = my_sphere.lats
    lons = my_sphere.lons

    def da_convert(arr):
        return xr.DataArray(arr, dims=['lat', 'lon'], coords={'lat':lats, 'lon':lons})

    '''
    ForcingFn takes one input which is the previous state. In this case there is no previous 
    state so we put zero
    '''

    forcing_real = my_sphere.spectogrd(perturb_actual)

    forcing_ds = da_convert(forcing_real)
    forcing_ds = forcing_ds * np.exp(-(forcing_ds.lat - lat_center)**2/dlat**2)
    forcing_ds = forcing_ds.values
    print(h['g'].shape, forcing_ds.shape)
    h['g'] = forcing_ds.T *meter

    print(abs(h['g']).max())

    #initial condition: saturation field
    delta_latQ = np.pi/3
    q_sat['g'] = Q * np.exp(-(lat_field['g'] - lat_center)**2 /(delta_latQ**2)) 


    # Create a trivial IVP containing the fields
    problem = d3.IVP([h, q_sat], namespace=locals())

    problem.add_equation("dt(h) = 0")
    problem.add_equation("dt(q_sat) = 0")

    # Create solver
    solver = problem.build_solver(d3.RK222)
    solver.stop_sim_time = 1
    # Checkpoint handler
    checkpoints = solver.evaluator.add_file_handler(
        data_path + f'IVP/{args.data_folder}/',
        sim_dt=1,
        max_writes=1
    )

    checkpoints.add_tasks(solver.state)

    try:
        logger.info('Starting main loop')
        while solver.proceed:
            solver.step(1)
            if (solver.iteration-1) % 10 == 0:
                logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, 1))
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        solver.log_stats()



