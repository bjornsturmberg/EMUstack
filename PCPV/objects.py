# doc
import os
import materials

data_location = '../PCPV/Data/'

class SolarCell(object):
    """ Represents Solar Cell structure
    """
    def __init__(self, period, radius1, radius2, radius3, radius4, ff, 
        mesh_file = 'NEED_FILE.geo',
        set_ff = False, height_1 = 2330, height_2 = 2330, num_h = 1,
        inclusion_a = materials.Si_c, inclusion_b = materials.Si_c, background = materials.Air,
        superstrate = materials.Air, substrate = materials.SiO2_a,
        loss = True, lx = 1, ly = 1, mesh_format = 1, nb_typ_el = 4, make_mesh_now = False,
        lc = 0.09, lc1 = 1.7, lc2 = 1.7, lc3 = 2.7, lc4 = 1.9, lc5 = 1.9, lc6 = 1.9):
        self.radius1     = radius1
        self.radius2     = radius2
        self.radius3     = radius3
        self.radius4     = radius4
        self.period      = period
        self.ff          = ff
        self.height_1    = height_1
        self.height_2    = height_2
        self.num_h       = num_h
        self.inclusion_a = inclusion_a
        self.inclusion_b = inclusion_b
        self.background  = background
        self.superstrate = superstrate
        self.substrate   = substrate
        self.loss        = loss
        self.lx          = lx
        self.ly          = ly
        self.nb_typ_el   = nb_typ_el
        self.set_ff      = set_ff
        if substrate == superstrate:
            self.has_substrate = 0
        else:
            self.has_substrate = 1
        if make_mesh_now:
            self.lc  = lc
            self.lc1 = lc1
            self.lc2 = lc2
            self.lc3 = lc3
            self.lc4 = lc4
            self.lc5 = lc5
            self.lc6 = lc6
            self.mesh_format = mesh_format
            self.make_mesh()
        else:
            self.mesh_file   = mesh_file
            self.mesh_format = mesh_format


    def make_mesh(self):

        if self.radius4 > 0:
            supercell = 4
            msh_name  =  '%(d)i_%(radius)i_%(radiuss)i_%(radiusss)i_%(radiussss)i' % {
           'd' : self.period, 'radius' : self.radius1, 'radiuss' : self.radius2, 
           'radiusss' : self.radius3,  'radiussss' : self.radius4}
        elif self.radius3 > 0:
            supercell = 3
            msh_name  =  '%(d)i_%(radius)i_%(radiuss)i_%(radiusss)i' % {
           'd' : self.period, 'radius' : self.radius1, 'radiuss' : self.radius2, 
           'radiusss' : self.radius3}
        elif self.radius2 > 0:
            supercell = 2
            if self.set_ff == False:
                msh_name  =  '%(d)i_%(radius)i_%(radiuss)i' % {
               'd' : self.period, 'radius' : self.radius1, 'radiuss' : self.radius2}
            else:
                msh_name  =  '%(d)i_%(radius)i_f_%(ff)i' % {
               'd' : self.period, 'radius' : self.radius1, 'ff' : 100*round(self.ff,2)}
        elif self.radius1 > 0:
            supercell = 1
            if self.set_ff == False:
                msh_name  =  '%(d)i_%(radius)i' % {'d' : self.period, 'radius' : self.radius1}
            else:
                msh_name  =  '%(d)i_f_%(ff)i' % {'d' : self.period, 'ff' : 100*round(self.ff,2)}
        else:
            # except KeyError:
            raise  NotImplementedError, "must have at least one cylinder of nonzero radius."

        self.mesh_file = msh_name + '.mail'

        if not os.path.exists(data_location + msh_name + '.geo'):
            geo_tmp = open(data_location + '%s_msh_template.geo' % supercell, "r").read()
            geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
            geo = geo.replace('d_in_nm = 0;', "d_in_nm = %i;" % self.period)
            geo = geo.replace('a1 = 0;', "a1 = %i;" % self.radius1)
            geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
            geo = geo.replace('lc1 = lc/1;', "lc1 = lc/%f;" % self.lc1)
            geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
            geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
            if supercell == 1 & self.set_ff == True:
                geo = geo.replace('radius1 = (a1/d_in_nm)*d;','radius1 = ((ff*(d)^2)/3.14159265)^0.5;')
            if supercell == 2 & self.set_ff == False:
                geo = geo.replace('radius2 = ((ff*(d)^2)/3.145 - (radius1^2))^0.5;', "radius2 = (%i/d_in_nm)*d;" % self.radius1)    
            if supercell > 1:
                geo = geo.replace('a2 = 0;', "a2 = %i;" % self.radius2)
                geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
            if supercell > 2:
                geo = geo.replace('a3 = 0;', "a3 = %i;" % self.radius3)
                geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)
            if supercell > 3:
                geo = geo.replace('a4 = 0;', "a4 = %i;" % self.radius4)
                geo = geo.replace('lc6 = lc/1;', "lc6 = lc/%f;" % self.lc6)

            open(data_location + msh_name + '.geo', "w").write(geo)
            
            gmsh_cmd = './'+ data_location + 'gmsh_conversion/' + "conv_gmsh_py.exe "+ data_location + msh_name
            os.system(gmsh_cmd)

     

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
    """
    0  debug         - 0 for full simulations, 1 for full details
    1  parallel_comp - Parallel computation 1, else single core
    1  Loss          - Include loss 1, else set Im(n) to zero
    1  hole_wire     - Choose holes 0, nanowires 1, uniform slab 2 or const 3
    0  TEmode        - Calculate t,r,a for the incident TE polarisation
    1  TMmode        - Calculate t,r,a for the incident TM polarisation
    1  traLambda     - Save t,r,a for each wavelength
    1  PropModes     - Save # of propagating modes & det(I-RPRP) 1, save SVD 2
    0  PrintSolution - Save FEM Bloch modes for lambda selected below.
    0  PrintSupModes - Save supermode field distributions
    1  PrintOmega    - Save dispersion data (beta for each lambda)
    0  PrintAll      - Save J overlap, orthogonal
    0  Checks        - Check completeness, energy conservation

    1                    # Selected formulation (1=E-Field, 2=H-Field)
    1                    # Number of wavelength to be scanned (n_lambda)
    315.0d0  699.0d0     # Wavelength range in nm (lambda_1, lambda_2)
    0                    # Sets wavelengths as in data files, (Lambda_Res)
    600                  # Lattice vector in nanometers (d_in_nm)
    2                    # Boundary conditions (0=Dirichlet,1=Neumann,2=Periodic)
    1.0d0                # X dimension of unit cell (lx)
    1.0d0                # Y dimension of unit cell (ly)
    350                  # Number of FEM Eigenvalues to be calculated (nval)
    725                  # Dimensions of reduced matrix (must be > 2*nval)(nvect)
    30                   # Maximum number of iterations for convergence (itermax)
    0.0d0                # ARPACK accuracy (0.0 for machine precision)(tol) 
    5.0d0 0.0d0          # Re and Im of parameter for shift-and-invert method
    """
    def __init__(self, debug = 0, traLambda = 1 , PrintOmega = 1, PrintSolution = 0, 
        PrintSupModes = 0, PrintAll = 0, Checks = 0, PropModes = 0, q_average = 0, 
        plot_real = 1, plot_imag = 0, plot_abs = 0, tol = 0, E_H_field = 1, 
        i_cond = 2, itermax = 30, incident = 0, what4incident = 2, out4incident = 0,
        Animate = False):
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
        self.Animate       = Animate
