import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

import dedalus_grid_setup

dtype = np.float64
#################################################################
# Simulation units

def gaussian_height_qsat_IC(H,Q, dlat, dlon, lat_center, lon_center):
    """
    H: resting depth in meters
    Q: max qsat in meters
    dlat: meridional extent of gaussian in degrees
    dlon: zonal extent of gaussian in degrees
    lat_center: gaussian center latitude in degrees
    lon_center: gaussian center longitude in degrees
    """
    ug = dedalus_grid_setup.units_and_grid(H, Q, 5000, 0.1) #units and grids #chi and tau fields are just placeholders
    #, not relevant to this module

    H = ug.H
    Q = ug.Q #scale of saturation vapor pressure in meters of water in the column (50mm water about how much is found in the equatorial atm)

  
    dist = ug.dist
    basis = ug.basis

    h = ug.h
    q_sat = ug.q_sat

    phi_field = ug.phi_field
    lat_field = ug.lat_field

    print('Setup done')
    print(f'beginning moist simulation')

    ##########################################################################

    # Initial conditions: gaussian bump
    lat_center, lon_center = np.deg2rad(lat_center), np.deg2rad(lon_center) #edit

    #disturbance extents
    delta_lat_dist = np.deg2rad(dlat)
    delta_lon_dist = np.deg2rad(dlon)

    #sat vap falloff
    delta_latQ = np.pi/3 #meridional gradient


    bump = dist.Field(name='bump', bases = basis) 
    bump['g'] = 1 * (np.exp(-(lat_field['g'] - lat_center)**2 /(delta_lat_dist**2)) *
                    np.exp(-(phi_field['g'] - lon_center)**2 /(delta_lon_dist**2)))


    #actually initial conditions: create gaussian geopotential bump, then BVP for psi 
    h['g'] = - H/30 * bump['g']  #just the perturbation, low pressure

    #initial condition: saturation field
    q_sat['g'] = Q * np.exp(-(lat_field['g'] - lat_center)**2 /(delta_latQ**2)) 

    return h['g'], q_sat['g']

    

