import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

from dedalus.core.operators import UnaryGridFunction
import dedalus_grid_setup

import argparse

data_path = '/home/rajat/Desktop/Project/data/'
dtype = np.float64
####################################################
parser = argparse.ArgumentParser(description="target save folder and spinup IC choice")

parser.add_argument("--run_days", type=int, required=True,
                    help="No. of simulation days to run in days, default = 15 days")
parser.add_argument("--spinup", type=int, required=True, 
                    help="start from spunup up state  or everyting zero except the IC. Should be 0(no spinup) or 1(yes spinup)")
parser.add_argument("--data_folder", type=str,required=True,
                    help="folder to which to save all data. Will also take the initial condition from the checkpoint files in this folder")

parser.add_argument("--IC_type", type=str, required=True,
                    help="random or gaussian initial condition. dlon and lon_center are useless for random option, just put anything")

parser.add_argument("--dlat", type=float, required=True, help="lat extent of gaussian in degrees")
parser.add_argument("--dlon", type=float, required=True, help="lon extent of gaussian in degrees")
parser.add_argument("--lat_center", type=float, default=0, help="lat center of gaussian in degrees")
parser.add_argument("--lon_center", type=float, default=180, help="lon center of gaussian in degrees")

parser.add_argument("--conv_time", type=float, default=2.4, help="convective adjustment timescale in hours")

args = parser.parse_args()

if (args.spinup != 0) & (args.spinup != 1): 
    raise ValueError("spinup input should be 0(no spinup) or 1(yes spinup)")

#################################################################
# Simulation units

convective_timescale = args.conv_time #in hours
convective_timescale = convective_timescale/24 #in days
ug = dedalus_grid_setup.units_and_grid(300, 0.05, 5000, convective_timescale) #units and grids
meter = ug.meter
hour = ug.hour
second = ug.second

Nphi = ug.Nphi  #2pi/512 = 0.7 degrees, 78km at equator, scales as sin(theta) # theta = colatitude = 90 - latitude
Ntheta = ug.Ntheta #pi/256 = 0.7 degrees, 78km
dealias = ug.dealias
R = ug.R #so Radius of earth = 1 unit. Meter is a fraction of that
Omega = ug.Omega
g = ug.g
H = ug.H

timestep = 450 * second
stop_sim_time = 24*args.run_days* hour

Q = ug.Q #scale of saturation vapor pressure in meters of water in the column (50mm water about how much is found in the equatorial atm)
chi = ug.chi #non dimensional conversion factor between q and h
tau = ug.tau #condensation/evaporation timescale

delx = 2*np.pi * R /Nphi
nu4 =1*delx**4
hypervis = lambda A: nu4*d3.lap(d3.lap(A))
zcross = lambda A: d3.MulCosine(d3.skew(A)) #multiply by cos theta, i.e. sin(lat). 
                                                #skew rotates a vector by 90 degrees, so (a,b) --> (-b, a)
    #so this is \hat(z) \cross \vec(A)


class DedalusHeaviside(UnaryGridFunction):#input shout be (q-q_sat or q_sat-q depending on P or E)
    """A wrapper to safely apply np.heaviside to grid-space data in Dedalus."""
    def __init__(self, arg, **kwargs):
        # We define a localized wrapper function compatible with NumPy ufuncs
        func = lambda x, out: np.heaviside(x, 0, out=out)
        super().__init__(func, arg, **kwargs)

############################################################################
coords = ug.coords
dist = ug.dist
basis = ug.basis

u = ug.u
ex = ug.ex
ey = ug.ey
h = ug.h
psi = ug.psi
kappa = ug.kappa
q = ug.q
P = ug.P
E = ug.E
q_sat = ug.q_sat
##############################################################################

if args.IC_type == 'gaussian':
    from gaussian_height_qsat_IC import gaussian_height_qsat_IC
    h['g'], q_sat['g'] = gaussian_height_qsat_IC(H/meter, Q/meter, args.dlat, args.dlon, 
                                                 args.lat_center, args.lon_center)
    print(f'{args.IC_type} set')


elif args.IC_type == 'random':
    try:
        problem = d3.IVP([h, q_sat], namespace=locals())
        problem.add_equation("dt(h) = 0")
        problem.add_equation("dt(q_sat) = 0")

        # Create solver
        solver = problem.build_solver(d3.RK222)
        solver.load_state(data_path + 
            f'IVP/{args.data_folder}/{args.data_folder}_s1.h5', index=-1,allow_missing=True)

        h['g'] = solver.state[0]['g']
        q_sat['g'] = solver.state[1]['g']

        print(f'{args.IC_type} set')
        
    except:
        raise Exception
        
else:
    raise ValueError('IC_type only permits "gaussian" or "random" ')


print(abs(h['g'].max())/meter, q_sat['g'].max()/meter)
q_sat.change_scales(1)
h.change_scales(1)

##########################################################################

print(f"Starting {args.run_days} day run")

# Problem
# -div(skew(A)) = k dot curl(A)
c = dist.Field(name='c')
d = dist.Field(name='d')
# Problem

problem = d3.IVP([u, psi, kappa, h, q, P, E, c, d], namespace=locals())
problem.add_equation("dt(u) + hypervis(u) = - u@grad(u) - g*grad(h) - 2*Omega*zcross(u)")
problem.add_equation("lap(psi) + c = -div(skew(u))")
problem.add_equation("lap(kappa) + d = div(u)") 
problem.add_equation("ave(psi) = 0")
problem.add_equation("ave(kappa) = 0")

problem.add_equation("dt(h) + hypervis(h) + H*div(u) = - div(h*u) - chi*(P-E)")

problem.add_equation("dt(q) + hypervis(q) = -P+E -  div(q*u) ")
problem.add_equation("P = (q-q_sat)*DedalusHeaviside(q - q_sat)/tau")
problem.add_equation("E = (q_sat-q)*DedalusHeaviside(q_sat - q)/tau")


# Solver
solver = problem.build_solver(d3.RK443)
solver.stop_sim_time = stop_sim_time

##############################################################################
print('Assigning initial state')
if args.spinup == 0:
    u['g']= 0
    psi['g'] = 0
    kappa['g'] = 0
    q['g'] = q_sat['g']
    P['g'] = 0
    E['g'] = 0

    foldername = args.data_folder + '/no_spinup'

if args.spinup == 1:
    solver.load_state(data_path + 
    f'IVP/{args.data_folder}/spinup_checkpoint/dry_part/dry_part_s1.h5', index=-1,allow_missing=True)

    solver.load_state(data_path + 
    f'IVP/{args.data_folder}/spinup_checkpoint/moist_part/moist_part_s1.h5', index=-1,allow_missing=True)

    foldername = args.data_folder
    
print(foldername)

##############################################################################

# Analysis
snapshots = solver.evaluator.add_file_handler(data_path + f'IVP/{foldername}/tau_c_{args.conv_time}_full_run', 
                                              sim_dt=6*hour, max_writes=16*5+1)#10days 3 hourly in one file

snapshots.add_task(d3.dot(u, ex)/(meter/second), name='u')
snapshots.add_task(-d3.dot(u, ey)/(meter/second), name='v') #note here that v as calculated is d\theta/dt, ie. colatitude velocity. so the actual meridional velocity is opposite that, positive northwards
#
snapshots.add_task(psi/(meter**2/second), name='psi')
snapshots.add_task(kappa/(meter**2/second), name='velpotn')
snapshots.add_task(-d3.div(d3.skew(u))/(1/second), name='vort')
snapshots.add_task(d3.div(u)/(1/second), name='divg')
snapshots.add_task(q/(meter), name='q')
snapshots.add_task(h/(meter), name='h')
snapshots.add_task(P/(meter/second), name='P')
snapshots.add_task(E/(meter/second), name='E')

# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        solver.step(timestep)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()

print('sim done')
