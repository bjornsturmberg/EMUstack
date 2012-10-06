# doc
import os
import numpy as np
import random
import materials
from calculate_ff     import calculate_ff

data_location = '../PCPV/Data/'

class NanoStruct(object):
    """ Represents structured layer
    """
    def __init__(self, period, ff, radius1, radius2=0, radius3=0, radius4=0, radius5=0,
        radius6=0, radius7=0, radius8=0, radius9=0, radius10=0, radius11=0, radius12=0, radius13=0,
        radius14=0, radius15=0, radius16=0, ellipticity = 0.0, square = False, mesh_file = 'NEED_FILE.geo',
        set_ff = False, ff_rand = False, height_1 = 2330, height_2 = 2330, num_h = 1,
        inclusion_a = materials.Si_c, inclusion_b = materials.Air, background = materials.Air,
        superstrate = materials.Air, substrate = materials.SiO2_a,
        loss = True, lx = 1, ly = 1, mesh_format = 1, nb_typ_el = 4, make_mesh_now = False, force_mesh = False,
        lc_bkg = 0.09, lc2_surf= 1.7, lc3_inc1 = 1.9, lc4_inc2 = 1.9, lc5_inc3 = 1.9, lc6_inc4 = 1.9,
        posx = 0, posy = 0, max_num_BMs = 100, label_nu = 0):
        self.radius1     = radius1
        self.radius2     = radius2
        self.radius3     = radius3
        self.radius4     = radius4
        self.radius5     = radius5
        self.radius6     = radius6
        self.radius7     = radius7
        self.radius8     = radius8
        self.radius9     = radius9
        self.radius10    = radius10
        self.radius11    = radius11
        self.radius12    = radius12
        self.radius13    = radius13
        self.radius14    = radius14
        self.radius15    = radius15
        self.radius16    = radius16
        self.period      = period
        self.ff          = ff
        self.ff_rand     = ff_rand
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
        self.max_num_BMs = max_num_BMs
        self.label_nu       = label_nu
        if substrate == superstrate:
            self.has_substrate = 0
        else:
            self.has_substrate = 1
        if make_mesh_now:
            self.lc  = lc_bkg
            self.lc2 = lc2_surf
            self.lc3 = lc3_inc1
            self.lc4 = lc4_inc2
            self.lc5 = lc5_inc3
            self.lc6 = lc6_inc4
            self.ellipticity = ellipticity
            if ellipticity > 1.0:
                raise TypeError, "ellipticity must be less than 1.0"
            self.square = square
            self.square_int = 0
            self.posx = posx
            self.posy = posy
            self.mesh_format = mesh_format
            self.force_mesh  = force_mesh
            self.make_mesh()
        else:
            self.mesh_file   = mesh_file
            self.mesh_format = mesh_format



    def make_mesh(self):

        if self.radius10 > 0:
            supercell = 16
            msh_name  =  '%(d)i_%(radius)i_%(radiuss)i_%(radiusss)i_%(radiussss)i_%(adiussss)i' % {
           'd' : self.period, 'radius' : self.radius1, 'radiuss' : self.radius2, 
           'radiusss' : self.radius3,'radiussss' : self.radius4, 'adiussss' : self.radius5}
        elif self.radius5 > 0:
            supercell = 9
            msh_name  =  '%(d)i_%(radius)i_%(radiuss)i_%(radiusss)i_%(radiussss)i_%(adiussss)i' % {
           'd' : self.period, 'radius' : self.radius1, 'radiuss' : self.radius2, 
           'radiusss' : self.radius3,'radiussss' : self.radius4, 'adiussss' : self.radius5}
        elif self.radius4 > 0:
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
                self.radius2 = ((self.ff*(self.period)**2)/3.14159265 - (self.radius1**2))**0.5
        elif self.radius1 > 0:
            supercell = 1
            if self.set_ff == False:
                msh_name  =  '%(d)i_%(radius)i' % {'d' : self.period, 'radius' : self.radius1}
            else:
                msh_name  =  '%(d)i_f_%(ff)i' % {'d' : self.period, 'ff' : 100*round(self.ff,2)}
                self.radius1 = np.sqrt((self.ff*(self.period)**2)/3.14159265)
        else:
            # except KeyError:
            raise  NotImplementedError, "must have at least one cylinder of nonzero radius."
 

        if self.ellipticity != 0:
            msh_name = msh_name + '_e_%(e)i' % {'e' : self.ellipticity*100,}
        if self.square == True:
            self.square_int = 1
            msh_name = msh_name + '_sq'
        if self.posx != 0:
            msh_name = msh_name + 'x%(e)i' % {'e' : self.posx*100,}
        if self.posy != 0:
            msh_name = msh_name + 'y%(e)i' % {'e' : self.posy*100,}

        self.mesh_file = msh_name + '.mail'
            
        # for blah in range(1,101,1):
        #     print blah
        #     msh_name = 'random_u_%i' % blah
        #     self.mesh_file = msh_name + '.mail'

        # msh_name = 'design-last_17'
        # self.mesh_file = msh_name + '.mail'

        if self.ff_rand == True:
            ff_tol = 0.0001
            min_a  = 50
            max_a  = (self.period/2.05)/np.sqrt(supercell)
            unit_period = (self.period/np.sqrt(supercell))
            mean = np.sqrt((self.ff*(unit_period)**2)/3.14159265)
            test_ff = 0
            while abs(test_ff-self.ff) > ff_tol:
                rad_array = []
                for i in range(supercell):
                    # stand_dev = 30
                    # select_radius = random.gauss(mean,stand_dev)
                    select_radius = random.uniform(min_a,max_a)
                    rad_array = np.append(rad_array,select_radius)

                test_ff = calculate_ff(self.period,rad_array[0],rad_array[1],rad_array[2],rad_array[3],rad_array[4],
                rad_array[5],rad_array[6],rad_array[7],rad_array[8],rad_array[9],rad_array[10],
                rad_array[11],rad_array[12],rad_array[13],rad_array[14],rad_array[15])
                print test_ff
                if supercell>3:
                    self.radius1   = rad_array[0]
                    self.radius2   = rad_array[1]
                    self.radius3   = rad_array[2]
                    self.radius4   = rad_array[3]
                if supercell>4:
                    self.radius5   = rad_array[4]
                    self.radius6   = rad_array[5]
                    self.radius7   = rad_array[6]
                    self.radius8   = rad_array[7]
                    self.radius9   = rad_array[8]
                if supercell>9:
                    self.radius10  = rad_array[9]
                    self.radius11  = rad_array[10]
                    self.radius12  = rad_array[11]
                    self.radius13  = rad_array[12]
                    self.radius14  = rad_array[13]
                    self.radius15  = rad_array[14]
                    self.radius16  = rad_array[15]
                test_ff = calculate_ff(self.period,rad_array[0],rad_array[1],rad_array[2],rad_array[3],rad_array[4],
                rad_array[5],rad_array[6],rad_array[7],rad_array[8],rad_array[9],rad_array[10],
                rad_array[11],rad_array[12],rad_array[13],rad_array[14],rad_array[15])

        if not os.path.exists(data_location + msh_name + '.mail') or self.force_mesh == True:
            geo_tmp = open(data_location + '%s_msh_template.geo' % supercell, "r").read()
            geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
            geo = geo.replace('d_in_nm = 0;', "d_in_nm = %i;" % self.period)
            geo = geo.replace('a1 = 0;', "a1 = %i;" % self.radius1)
            geo = geo.replace('ellipticity = 0;', "ellipticity = %f;" % self.ellipticity)
            geo = geo.replace('square = 0;', "square = %i;" % self.square)
            geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
            geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
            geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
            if self.posx != 0:
                geo = geo.replace('posx = 0;', "posx = %f;" % (self.posx/self.period*(self.period/(2*np.sqrt(supercell)) - self.radius1)))
            if self.posy != 0:
                geo = geo.replace('posy = 0;', "posy = %f;" % (self.posy/self.period*(self.period/(2*np.sqrt(supercell)) - self.radius1)))
            if supercell == 2:
                if self.set_ff == False:
                    geo = geo.replace('radius2 = ((ff*(d)^2)/3.14159265 - (radius1^2))^0.5;', "radius2 = (%i/d_in_nm)*d;" % self.radius2)    
            if supercell > 1:
                geo = geo.replace('a2 = 0;', "a2 = %i;" % self.radius2)
                geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
            if supercell > 2:
                geo = geo.replace('a3 = 0;', "a3 = %i;" % self.radius3)
                geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)
            if supercell > 3:
                geo = geo.replace('a4 = 0;', "a4 = %i;" % self.radius4)
                geo = geo.replace('lc6 = lc/1;', "lc6 = lc/%f;" % self.lc6)
            if supercell > 4:
                geo = geo.replace('a5 = 0;', "a5 = %i;" % self.radius5)
                geo = geo.replace('a6 = 0;', "a6 = %i;" % self.radius6)
                geo = geo.replace('a7 = 0;', "a7 = %i;" % self.radius7)
                geo = geo.replace('a8 = 0;', "a8 = %i;" % self.radius8)
                geo = geo.replace('a9 = 0;', "a9 = %i;" % self.radius9)
            if supercell > 9:
                geo = geo.replace('a10 = 0;', "a10 = %i;" % self.radius10)
                geo = geo.replace('a11 = 0;', "a11 = %i;" % self.radius11)
                geo = geo.replace('a12 = 0;', "a12 = %i;" % self.radius12)
                geo = geo.replace('a13 = 0;', "a13 = %i;" % self.radius13)
                geo = geo.replace('a14 = 0;', "a14 = %i;" % self.radius14)
                geo = geo.replace('a15 = 0;', "a15 = %i;" % self.radius15)
                geo = geo.replace('a16 = 0;', "a16 = %i;" % self.radius16)


            open(data_location + msh_name + '.geo', "w").write(geo)
            
            gmsh_cmd = './'+ data_location + 'gmsh_conversion/' + "conv_gmsh_py.exe "+ data_location + msh_name
            os.system(gmsh_cmd)
            # gmsh_cmd = 'gmsh '+ data_location + msh_name + '.msh'
            # gmsh_cmd = 'gmsh '+ data_location + msh_name + '.geo'
            # os.system(gmsh_cmd)
         
class ThinFilm(object):
    """ Represents homogeneous film """
    def __init__(self, simo_period, height_1 = 2330, height_2 = 2330,
         num_h = 1, film_material = materials.Si_c, 
        superstrate = materials.Air, substrate = materials.Air, 
        nu_tot_ords = 0, nu_prop_ords = 0, zero_ord = 0, 
        set_ord_in = 0, set_ord_out = 0, loss = True, label_nu = 0):
        self.simo_period   = simo_period
        self.height_1      = height_1
        self.height_2      = height_2
        self.num_h         = num_h
        self.film_material = film_material
        self.superstrate   = superstrate
        self.substrate     = substrate
        self.nu_tot_ords   = nu_tot_ords
        self.nu_prop_ords  = nu_prop_ords
        self.zero_ord      = zero_ord
        self.set_ord_in    = set_ord_in
        self.set_ord_out   = set_ord_out
        self.loss          = loss
        self.label_nu      = label_nu
       


class Light(object):
    """ Represents the light incident on structure """
    def __init__(self, Lambda, pol = 'TM', theta = 0, phi = 0):
        self.Lambda = Lambda
        if  pol == 'Unpol':
            self.pol = 0
        elif  pol == 'TE':
            self.pol = 1
        elif pol == 'TM':
            self.pol = 2
        elif pol == 'Left':
            self.pol = 3
        elif pol == 'Right':
            self.pol = 4
        elif pol == 'CD':
            self.pol = 5
        else:
            raise TypeError, "Polarisation must be either 'Unpol','TE', 'TM', 'Left', 'Right', or 'CD'."
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
        i_cond = 2, itermax = 30, incident = 0, 
        x_order_in = 0, x_order_out = 0, y_order_in = 0, y_order_out = 0,
        what4incident = 2, out4incident = 0, Animate = False, 
        var_BM_min = 370, max_order_PWs = 3, num_cores = 8, leave_cpus = False):
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
        self.x_order_in    = x_order_in
        self.x_order_out   = x_order_out
        self.y_order_in    = y_order_in
        self.y_order_out   = y_order_out
        self.what4incident = what4incident
        self.out4incident  = out4incident
        self.Animate       = Animate
        self.var_BM_min    = var_BM_min
        self.max_order_PWs = max_order_PWs
        self.num_cores     = num_cores
        if leave_cpus == True:
            # number of cpus to leave free
            import multiprocessing   as mp
            num_cores_sahand_gets = 8
            self.num_cores = mp.cpu_count() - num_cores_sahand_gets
