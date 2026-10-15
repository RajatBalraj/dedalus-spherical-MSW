import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

dtype = np.float64

class units_and_grid(object):
    def __init__(self, H, Q, chi, tau): #in meters, meters, dimensionless and day units respectively

        self.meter = 1 / 6.37122e6
        self.hour = 1
        self.second = self.hour/3600
        self.Nphi = 512
        self.Ntheta = 256
        self.dealias = 3/2
        self.R = 6.37122e6 * self.meter
        self.Omega = 7.292e-5 / self.second
        self.g = 9.80616 * self.meter / self.second**2

        self.H = H*self.meter
        self.Q = Q*self.meter
        self.chi = chi
        self.tau = tau*24*self.hour

        # Bases
        coords = d3.S2Coordinates('phi', 'theta') #azimuth and colatitude(90 - latitude)
        dist = d3.Distributor(coords, dtype=dtype)
        basis = d3.SphereBasis(coords, (self.Nphi, self.Ntheta), radius=self.R, dealias=self.dealias, dtype=dtype)


        # Fields
        u = dist.VectorField(coords, name='u', bases=basis) # u['g'] is a vector field has the shape 
        #(2 (two vector components u and v), 384(no. of orders(m, lon)), 192(no. of degrees of harmonics considered(l), lat)), which correspond to

        ex = dist.VectorField(coords, name='ex', bases=basis)
        ex['g'][0] = 1
        ex['g'][1] = 0

        ey = dist.VectorField(coords, name='ey', bases=basis)
        ey['g'][0] = 0
        ey['g'][1] = 1
        #this is because to find velocity field from height, we need two equations, we cant use just the laplacian height/divergence equation
        #either we use both the momentum eqns or the divergence as well as the vorticity eqations. I suppose this method is better.

        h = dist.Field(name='h', bases=basis) #h['g'] is a scalar field, shape (256, 128)

        psi = dist.Field(name='psi', bases=basis) #streamfunction
        kappa = dist.Field(name='kappa', bases=basis) #velocity potential

        q = dist.Field(name = 'q', bases = basis) #in units of m of water in atm column. Not appropriate to convert to pressure
        P = dist.Field(name='P', bases=basis) #precipitation
        E = dist.Field(name='E', bases=basis) #evaportion

        q_sat = dist.Field(name='q_sat', bases=basis) #time constant prescribed field

        phi, theta = dist.local_grids(basis)

        theta_field = dist.Field(name='theta', bases = basis)
        theta_field['g'] = dist.local_grids(basis)[1] + 0*dist.local_grids(basis)[0] #for making initial condition

        f = dist.Field(name='f', bases = basis) 
        f['g'] = 2* self.Omega * np.cos(theta_field['g'])

        phi_field = dist.Field(name='phi', bases = basis)
        phi_field['g'] = 0*dist.local_grids(basis)[1] + dist.local_grids(basis)[0]

        lat = np.pi / 2 - theta + 0*phi #adding the phi increases the dimension of lat # because theta is colatitude

        lat_field = dist.Field(name='lat', bases = basis)
        lat_field['g'] = np.pi/2 - theta_field['g'] 

###########################################################################################
        self.coords = coords
        self.dist = dist
        self.basis = basis

        self.u = u
        self.ex = ex
        self.ey = ey
        self.h = h
        self.psi = psi
        self.kappa = kappa
        self.q = q
        self.P = P
        self.E = E
        self.q_sat = q_sat
        self.phi = phi
        self.theta = theta
        self.f = f
        self.phi_field = phi_field
        self.theta_field = theta_field
        self.lat = lat
        self.lat_field = lat_field
