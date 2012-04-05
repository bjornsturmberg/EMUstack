# doc

import materials

class SolarCell(object):
    """ Represents Solar Cell structure
    """
    def __init__(self, radius1, radius2, period, ff,
        mesh_file, height = 2330,
        inclusion = materials.cSi, background = materials.air,
        superstrate = materials.air, substrate = materials.SiO2,
        loss = True, lx = 1, ly = 1, mesh_format = 1,
        nb_typ_el = 4, make_mesh_now = False):
        self.radius1     = radius1
        self.radius2     = radius2
        self.period      = period
        self.ff          = ff
        self.height      = height
        self.inclusion   = inclusion
        self.background  = background
        self.superstrate = superstrate
        self.substrate   = substrate
        self.loss        = loss
        self.lx          = lx
        self.ly          = ly
        self.mesh_format = mesh_format
        self.nb_typ_el   = nb_typ_el
        if substrate == superstrate:
            self.has_substrate = 0
        else:
            self.has_substrate = 1
        if make_mesh_now:
            self.make_mesh()
        else:
            self.mesh_file = mesh_file


    def make_mesh(self, mesh_res,):

        #self.mesh_file = ???
        pass

     

class Light(object):
    """ Represents the light incident on a PC.
    """
    def __init__(self, Lambda, pol = 'TM', theta = 0, phi = 0):
        """See docstring for Light"""
        self.Lambda = Lambda
        if  pol == 'TE':
            self.pol = 1
        elif pol == 'TM':
            self.pol = 2
        else:
            raise TypeError, "Polarisation must be either 'TE' or 'TM'."
        self.theta = theta
        self.phi   = phi

        

class Controls(object):
    """docstring for Controls"""
    def __init__(self, debug = 0, traLambda = 1 , PrintOmega = 1, PrintSolution = 0, 
        PrintSupModes = 0, PrintAll = 0, Checks = 0, PropModes = 0, q_average = 0, 
        plot_real = 1, plot_imag = 0, plot_abs = 0, tol = 0, E_H_field = 1, 
        i_cond = 2, itermax = 30, incident = 0, what4incident = 2, out4incident = 0):
        self.debug         = debug
        self.traLambda     = traLambda
        self.PrintOmega    = PrintOmega
        self.PrintSolution = PrintSolution
        self.PrintSupModes = PrintSupModes
        self.PrintAll      = PrintAll
        self.Checks        = Checks
        self.PropModes     = PropModes
        self.q_average     = q_average
        self.plot_real     = plot_real
        self.plot_imag     = plot_imag
        self.plot_abs      = plot_abs
        self.tol           = tol
        self.E_H_field     = E_H_field
        self.i_cond        = i_cond
        self.itermax       = itermax
        self.incident      = incident
        self.what4incident = what4incident
        self.out4incident  = out4incident