"""
    objects.py is a subroutine of EMUstack that contains the NanoStruct,
    ThinFilm and Light objects. These represent the properties of a 
    structured layer, a homogeneous layer and the incident light
    respectively.

    Copyright (C) 2013  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

    EMUstack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import numpy as np
import random
import materials
from mode_calcs import Simmo, Anallo

data_location = '../backend/data/'

# Acknowledgements
print '\n#################################################################\n' + \
      'EMUstack is brought to you by Bjorn Sturmberg, Kokou Dossou, \n' + \
      'Felix Lawrence and Lindsay Botton, with support from CUDOS and ARENA\n' + \
      'Starting EMUstack calculation ...\n' + \
      '#################################################################\n'




class NanoStruct(object):
    """ Represents a structured layer.

        INPUTS:

        - 'geometry'      : Either 1D or 2D structure; '1D_array', '2D_array'.

        - 'period'        : The diameter the unit cell in nanometers.

        - 'diameter1'     : The diameter of the inclusion in nm.
        - 'diameter2-16'  : The diameter of further inclusions in nm.

        - 'height_nm'     : The thickness of the layer in nm or 'semi_inf'
            for a semi-infinte layer.

        - 'inclusion_a'   : A :Material: instance for first inclusion, 
            specified as dispersive refractive index (eg. materials.Si_c) 
            or nondispersive complex number (eg. Material(1.0 + 0.0j)).

        - 'inclusion_b'   : "  " for the second inclusion medium.
        - 'background'    : "  " for the background medium.

        - 'loss'          : If False, Im(n) = 0, if True n as in :Material: instance.

        - 'hyperbolic'    : If True FEM looks for Eigenvalues around n**2 * k_0**2
            rather than the regular n**2 * k_0**2 - alpha**2 - beta**2.

        - 'ellipticity'   : If != 0, inclusion has given ellipticity, with b=diameter,
           a=diameter-ellipticity*diameter. NOTE: only implemented for 1 inclusion.

        - 'square'        : If True, '2D_array' has square NWs (ie. 2D grating).

        - 'ff'            : The fill fraction of the inclusions. If non-zero, 
            the specified diameters are overritten s.t. given ff is achieved,
            otherwise ff is calculated from parameters and stored in self.ff.
        - 'ff_rand'       : If True, diameters overritten with random diameters,
            s.t. the ff is as assigned. Must provide non-zero dummy diameters.
        
        - 'posx'          : Shift NWs laterally towards center (each other), 
            posx is a fraction of the distance possible before NWs touching.
        - 'posy'          : Shift NWs vertically "  ". 
 
        - 'small_d'       : Distance between 2 inclusions of interleaved 1D grating.

        - 'make_mesh_now' : If True, program creates a FEM mesh with provided 
            NanoStruct parameters. If False, must provide mesh_file name of 
            existing .mail that will be run despite NanoStruct parameters.

        - 'force_mesh'    : If True, a new mesh is created despite existance of 
            mesh with same parameter. This is used to make mesh with equal 
            period etc. but different lc refinement.

        - 'mesh_file'     : If using a set premade mesh give its name including 
            .mail (eg. 600_60.mail), it must be located in backend/Data/

        - 'lc_bkg'        : Length constant of meshing of background medium.
            (smaller = finer mesh)
        - 'lc2'           : factor by which lc_bkg should be reduced on inclusion 
            surfaces; lc_surface = cl_bkg / lc2.
        - 'lc3-6'         : "  " at center of inclusions.

        - `plot_modes'    : Plot modes (ie. FEM solutions) in gmsh format, 
            you get epsilon*|E|^2 & either real/imag/abs of 
            x,y,z components, field vectors.

        - `plot_real'     : Plot the real part of modal fields.
        - `plot_imag'     : Plot the imaginary part of modal fields.
        - `plot_abs'      : Plot the absolute value of modal fields.
        - `plotting3d'    : Plot the fields in 3D.
    """


    def __init__(self, geometry, period, diameter1, 
        height_nm=2330,
        inclusion_a = materials.Si_c, inclusion_b = materials.Air,
        background = materials.Air, loss=True, hyperbolic=False,
        diameter2=0,  diameter3=0, diameter4=0, diameter5=0, diameter6=0, 
        diameter7=0, diameter8=0, diameter9=0, diameter10=0, diameter11=0, 
        diameter12=0, diameter13=0,diameter14=0, diameter15=0, diameter16=0, 
        ellipticity=0.0, square=False,
        ff=0, ff_rand=False, 
        small_d=0, posx=0, posy=0,
        make_mesh_now=True, force_mesh=False, 
        mesh_file='NEED_FILE.mail', 
        lc_bkg=0.09, lc2=1.0, lc3=1.0, lc4=1.0, lc5=1.0, lc6=1.0,
        plot_modes=False, plot_real=1, plot_imag=0, plot_abs=0,
        plotting3d = False):
        self.geometry      = geometry
        self.period        = period
        self.diameter1     = diameter1
        self.height_nm     = height_nm
        self.inclusion_a   = inclusion_a
        self.inclusion_b   = inclusion_b
        self.background    = background
        self.loss          = loss
        self.hyperbolic    = hyperbolic
        self.diameter2     = diameter2
        self.diameter3     = diameter3
        self.diameter4     = diameter4
        self.diameter5     = diameter5
        self.diameter6     = diameter6
        self.diameter7     = diameter7
        self.diameter8     = diameter8
        self.diameter9     = diameter9
        self.diameter10    = diameter10
        self.diameter11    = diameter11
        self.diameter12    = diameter12
        self.diameter13    = diameter13
        self.diameter14    = diameter14
        self.diameter15    = diameter15
        self.diameter16    = diameter16
        self.ellipticity   = ellipticity
        if ellipticity > 1.0:
            raise TypeError, "ellipticity must be less than 1.0"
        self.square        = square
        if square == True:
            self.square_int = 1
        else:
            self.square_int = 0
        if ff == 0:
            if geometry == '2D_array':
                self.ff = calculate_ff(square,period,diameter1,diameter2,
                    diameter3,diameter4,diameter5,diameter6,diameter7,diameter8,diameter9,
                    diameter10,diameter11,diameter12,diameter13,diameter14,diameter15,
                    diameter16,ellipticity)
            elif geometry == '1D_array':
                self.ff        = (diameter1 + diameter2)/period
        else:
            self.ff = ff
            if diameter2 != 0:
                self.diameter2 = 2*((ff*(period)**2)/np.pi - ((diameter1/2)**2))**0.5
            else:
                self.diameter1 = 2*np.sqrt((ff*(period)**2)/np.pi)
        self.ff_rand       = ff_rand
        self.small_d       = small_d
        self.posx          = posx
        self.posy          = posy
        self.force_mesh    = force_mesh
        self.lc            = lc_bkg
        self.lc2           = lc2
        self.lc3           = lc3
        self.lc4           = lc4
        self.lc5           = lc5
        self.lc6           = lc6
        if make_mesh_now == True:
            self.make_mesh()
        else:
            self.mesh_file = mesh_file
        if plot_modes == True:
            self.plot_modes = 1
            if not os.path.exists("Bloch_Fields"): os.mkdir("Bloch_Fields")
            if not os.path.exists("Bloch_Fields/PNG"): os.mkdir("Bloch_Fields/PNG")
        else: self.plot_modes = 0
        self.plot_real     = plot_real
        self.plot_imag     = plot_imag
        self.plot_abs      = plot_abs 
        self.plotting3d    = plotting3d 

    def make_mesh(self):
        if self.geometry == '2D_array':
            if self.diameter10 > 0:
                supercell = 16
                msh_name  =  '%(d)s_%(diameter)s_%(diameters)s_%(diameterss)s_%(diametersss)s_%(adiussss)s' % {
               'd' : dec_float_str(self.period), 'diameter' : dec_float_str(self.diameter1), 
               'diameters' : dec_float_str(self.diameter2), 'diameters' : dec_float_str(self.diameter2), 
               'diameterss' : dec_float_str(self.diameter3),'diametersss' : dec_float_str(self.diameter4), 
               'adiussss' : dec_float_str(self.diameter5)}
            elif self.diameter5 > 0:
                supercell = 9
                msh_name  =  '%(d)s_%(diameter)s_%(diameters)s_%(diameterss)s_%(diametersss)s_%(adiussss)s' % {
               'd' : dec_float_str(self.period), 'diameter' : dec_float_str(self.diameter1), 
               'diameters' : dec_float_str(self.diameter2), 'diameterss' : dec_float_str(self.diameter3),
               'diametersss' : dec_float_str(self.diameter4), 'adiussss' : dec_float_str(self.diameter5)}
            elif self.diameter4 > 0:
                supercell = 4
                msh_name  =  '%(d)s_%(diameter)s_%(diameters)s_%(diameterss)s_%(diametersss)s' % {
               'd' : dec_float_str(self.period), 'diameter' : dec_float_str(self.diameter1), 
               'diameters' : dec_float_str(self.diameter2), 'diameterss' : dec_float_str(self.diameter3), 
               'diametersss' : dec_float_str(self.diameter4)}
            elif self.diameter3 > 0:
                supercell = 3
                msh_name  =  '%(d)s_%(diameter)s_%(diameters)s_%(diameterss)s' % {
               'd' : dec_float_str(self.period), 'diameter' : dec_float_str(self.diameter1),
               'diameters' : dec_float_str(self.diameter2), 'diameterss' : dec_float_str(self.diameter3)}
            elif self.diameter2 > 0:
                supercell = 2
                msh_name  =  '%(d)s_%(diameter)s_%(diameters)s' % {'d' : dec_float_str(self.period), 
                'diameter' : dec_float_str(self.diameter1), 'diameters' : dec_float_str(self.diameter2)}
            elif self.diameter1 > 0:
                supercell = 1
                msh_name  =  '%(d)s_%(diameter)s' % {'d' : dec_float_str(self.period), 
                'diameter' : dec_float_str(self.diameter1)}
            else:
                raise ValueError, "must have at least one cylinder of nonzero diameter."

            if self.ellipticity != 0:
                msh_name = msh_name + '_e_%(e)s' % {'e' : dec_float_str(self.ellipticity),}
            if self.square == True:
                msh_name = msh_name + '_sq'
            if self.posx != 0:
                msh_name = msh_name + 'x%(e)s' % {'e' : dec_float_str(self.posx),}
            if self.posy != 0:
                msh_name = msh_name + 'y%(e)s' % {'e' : dec_float_str(self.posy),}

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
                max_a  = (self.period/1.05)/np.sqrt(supercell)
                unit_period = (self.period/np.sqrt(supercell))
                mean = np.sqrt((self.ff*(unit_period)**2)/np.pi)
                test_ff = 0
                while abs(test_ff-self.ff) > ff_tol:
                    rad_array = []
                    for i in range(supercell):
                        # stand_dev = 30
                        # select_diameter = random.gauss(mean,stand_dev)
                        select_diameter = random.uniform(min_a,max_a)
                        rad_array = np.append(rad_array,select_diameter)

                    test_ff = calculate_ff(self.square, self.period,rad_array[0],rad_array[1],rad_array[2],rad_array[3],rad_array[4],
                    rad_array[5],rad_array[6],rad_array[7],rad_array[8],rad_array[9],rad_array[10],
                    rad_array[11],rad_array[12],rad_array[13],rad_array[14],rad_array[15])
                    print test_ff
                    if supercell>3:
                        self.diameter1   = rad_array[0]
                        self.diameter2   = rad_array[1]
                        self.diameter3   = rad_array[2]
                        self.diameter4   = rad_array[3]
                    if supercell>4:
                        self.diameter5   = rad_array[4]
                        self.diameter6   = rad_array[5]
                        self.diameter7   = rad_array[6]
                        self.diameter8   = rad_array[7]
                        self.diameter9   = rad_array[8]
                    if supercell>9:
                        self.diameter10  = rad_array[9]
                        self.diameter11  = rad_array[10]
                        self.diameter12  = rad_array[11]
                        self.diameter13  = rad_array[12]
                        self.diameter14  = rad_array[13]
                        self.diameter15  = rad_array[14]
                        self.diameter16  = rad_array[15]
                    test_ff = calculate_ff(self.square, self.period,rad_array[0],rad_array[1],rad_array[2],rad_array[3],rad_array[4],
                    rad_array[5],rad_array[6],rad_array[7],rad_array[8],rad_array[9],rad_array[10],
                    rad_array[11],rad_array[12],rad_array[13],rad_array[14],rad_array[15])




            if not os.path.exists(data_location + msh_name + '.mail') or self.force_mesh == True:
                geo_tmp = open(data_location + '%s_msh_template.geo' % supercell, "r").read()
                geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                geo = geo.replace('d_in_nm = 0;', "d_in_nm = %f;" % self.period)
                geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                geo = geo.replace('ellipticity = 0;', "ellipticity = %f;" % self.ellipticity)
                geo = geo.replace('square = 0;', "square = %i;" % self.square_int)
                geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                if self.posx != 0:
                    # appropriate for old definition of fraction of distance to touching
                    geo = geo.replace('posx = 0;', "posx = %f;" % (self.posx/self.period*(self.period/(2*np.sqrt(supercell)) - self.diameter1/2.0)))
                    # appropriate for % shift of distance of centre point to (ind) unitcell boundary (ie d/2)
                    # geo = geo.replace('posx = 0;', "posx = %f;" % float(self.posx/supercell))
                if self.posy != 0:
                    geo = geo.replace('posy = 0;', "posy = %f;" % (self.posy/self.period*(self.period/(2*np.sqrt(supercell)) - self.diameter1/2.0)))
                    # geo = geo.replace('posy = 0;', "posy = %f;" % float(self.posy/supercell))
                if supercell > 1:
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                if supercell > 2:
                    geo = geo.replace('a3 = 0;', "a3 = %f;" % self.diameter3)
                    geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)
                if supercell > 3:
                    geo = geo.replace('a4 = 0;', "a4 = %f;" % self.diameter4)
                    geo = geo.replace('lc6 = lc/1;', "lc6 = lc/%f;" % self.lc6)
                if supercell > 4:
                    geo = geo.replace('a5 = 0;', "a5 = %f;" % self.diameter5)
                    geo = geo.replace('a6 = 0;', "a6 = %f;" % self.diameter6)
                    geo = geo.replace('a7 = 0;', "a7 = %f;" % self.diameter7)
                    geo = geo.replace('a8 = 0;', "a8 = %f;" % self.diameter8)
                    geo = geo.replace('a9 = 0;', "a9 = %f;" % self.diameter9)
                if supercell > 9:
                    geo = geo.replace('a10 = 0;', "a10 = %f;" % self.diameter10)
                    geo = geo.replace('a11 = 0;', "a11 = %f;" % self.diameter11)
                    geo = geo.replace('a12 = 0;', "a12 = %f;" % self.diameter12)
                    geo = geo.replace('a13 = 0;', "a13 = %f;" % self.diameter13)
                    geo = geo.replace('a14 = 0;', "a14 = %f;" % self.diameter14)
                    geo = geo.replace('a15 = 0;', "a15 = %f;" % self.diameter15)
                    geo = geo.replace('a16 = 0;', "a16 = %f;" % self.diameter16)

                open(data_location + msh_name + '.geo', "w").write(geo)
                
                gmsh_cmd = './'+data_location+'gmsh_conversion/'+"conv_gmsh_py.exe "+data_location+msh_name
                os.system(gmsh_cmd)

            # Automatically show created mesh in gmsh.
                # gmsh_cmd = 'gmsh '+ data_location + msh_name + '.msh'
                # os.system(gmsh_cmd)
                # gmsh_cmd = 'gmsh '+ data_location + msh_name + '.geo'
                # os.system(gmsh_cmd)


        elif self.geometry == '1D_array':
            if self.diameter2 > 0:
                supercell = 2
                msh_name  =  '1D_%(d)s_%(diameter)s_%(diameters)s' % {
               'd' : dec_float_str(self.period), 'diameter' : dec_float_str(self.diameter1),
               'diameters' : dec_float_str(self.diameter2)}
            elif self.diameter1 > 0:
                supercell = 1
                msh_name  =  '1D_%(d)s_%(diameter)s' % {'d' : dec_float_str(self.period), 
                    'diameter' : dec_float_str(self.diameter1)}
            else:
                raise ValueError, "must have at least one grating of nonzero width."

            self.mesh_file = msh_name + '.mail'    

            
            if not os.path.exists(data_location + msh_name + '.mail') or self.force_mesh == True:
                geo_tmp = open(data_location + '1D_%s_msh_template.geo' % supercell, "r").read()
                geo = geo_tmp.replace('d_in_nm = 0;', "d_in_nm = %f;" % self.period)
                geo = geo.replace('w1 = 0;', "w1 = %f;" % self.diameter1)
                geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                if supercell > 1:
                    geo = geo.replace('w2 = 0;', "w2 = %f;" % self.diameter2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                    geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                if self.small_d != 0:
                    # small distance between centre of gratings in nm
                    # calc complementary large distance, which is added to top & bottom
                    large_d_on_2 = (self.period - self.diameter1/2 - self.diameter2/2 - self.small_d)/2
                    posx1 = large_d_on_2 + self.diameter1/2
                    posx2 = large_d_on_2 + self.diameter2/2
                    posx3 = large_d_on_2 + self.diameter1 + ((self.small_d - self.diameter1/2 - self.diameter2/2)/2)
                    geo = geo.replace('posx1 = hy/4;', "posx1 = %f/d_in_nm;" % posx1)
                    geo = geo.replace('posx2 = hy/4;', "posx2 = %f/d_in_nm;" % posx2)
                    geo = geo.replace('posx3 = hy/2;', "posx3 = %f/d_in_nm;" % posx3)
                # if supercell > 1:
                #     geo = geo.replace('a2 = 0;', "a2 = %i;" % self.diameter2)
                #     geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                # if supercell > 2:
                #     geo = geo.replace('a3 = 0;', "a3 = %i;" % self.diameter3)
                #     geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)


                open(data_location + msh_name + '.geo', "w").write(geo)
                
                gmsh_cmd = './'+ data_location + 'gmsh_conversion/' + "conv_gmsh_py.exe "+ data_location + msh_name
                os.system(gmsh_cmd)
                # gmsh_cmd = 'gmsh '+ data_location + msh_name + '.msh'
                # gmsh_cmd = 'gmsh '+ data_location + msh_name + '.geo'
                # os.system(gmsh_cmd)
        else:
            raise ValueError, "Must be simulating either a '1D_array' or a '2D_array'."

    def calc_modes(self, light, **args):
        """ Run a simulation to find the structure's modes.

            INPUTS:

            - `light`: A :Light: instance representing incident light.

            - `args` : A dict of options to pass to :Simmo.calc_modes:.
        """
        simmo = Simmo(self, light) 

        simmo.calc_modes(**args)
        return simmo














class ThinFilm(object):
    """ Represents an unstructured homogeneous film.

        INPUTS:

        - 'period'         : Artificial period imposed on homogeneous film 
            to give consistently defined plane waves in terms of 
            diffraction orders of structured layers.

        - 'height_nm'      : The thickness of the layer in nm or 'semi_inf'
            for a semi-infinte layer.

        - 'material'       : A :Material: instance specifying the n of 
            the layer and related methods.

        - 'num_pw_per_pol' : Number of plane waves per polarisation.

        - 'loss'           : If False sets Im(n) = 0, if True leaves n as is.
    """
    def __init__(self, period, height_nm = 2330, material = materials.Si_c, 
        num_pw_per_pol = 0, loss = True):
        self.period         = period
        self.height_nm      = height_nm
        self.material       = material
        self.num_pw_per_pol = num_pw_per_pol
        self.loss           = loss
    
    def calc_modes(self, light):
        an = Anallo(self, light)
        an.calc_modes()
        return an


















class Light(object):
    """ Represents the light incident on structure.

        Incident angles may either be specified by `k_parallel` or by
        incident angles `theta` and `phi`, together with the refractive
        index `n_inc` of the incident medium.

        `wl_nm` and `k_pll` are both in unnormalised units.

        INPUTS:

        - `wl_nm`         : Wavelength, in nanometers.

        - `max_order_PWs` : Maximum plane wave order to include.

        - `k_parallel`    : The wave vector components (k_x, k_y)
            parallel to the interface planes. Units of nm^-1.

        - `theta`         : Polar angle of incidence in degrees.

        - `phi`           : Azimuthal angle of incidence in degrees.
    """
    def __init__(self, wl_nm, max_order_PWs = 2, k_parallel = [0.,0.], 
        theta = None, phi = None, n_inc = 1.):
        self.wl_nm = float(wl_nm)
        self._air_anallos = {}
        self.max_order_PWs  = max_order_PWs

        if None == theta and [0,0] == k_parallel:
            raise ValueError, "Specify incident angle either by \n\
            k_parallel OR by theta, phi and n_inc."

        if None == theta:
            self.k_pll = np.array(k_parallel, dtype='float64')
        else:
            # Check for inconsistent input
            if [0,0] != k_parallel or phi == None:
                raise ValueError, "Specify incident angle either by \n\
            k_parallel OR by theta, phi and n_inc."
            # Avoid the degeneracies that occur at normal incidence (FEM does not deal well with them)
            if abs(theta) < 1e-4: theta += 1e-4
            if abs(phi) < 1e-4: phi += 1e-4
            # Calculate k_parallel from incident angles
            k = 2 * np.pi * np.real(n_inc) / self.wl_nm
            theta *= np.pi / 180
            phi *= np.pi / 180
            self.k_pll = k*np.sin(theta) * np.array(
                        [np.cos(phi), np.sin(phi)], dtype='float64')


    def _air_ref(self, period):
        """ Return a :Anallo: corresponding to this :Light: in free space.

            The :Anallo: will have `len(anallo.k_z) == 2 * num_pw'.
        """

        if (period) in self._air_anallos:
            return self._air_anallos[(period)]
        else:
            air = ThinFilm(period = period, material = materials.Air)
            an = Anallo(air, self)

            an.is_air_ref = True

            kz = an.calc_kz()

            an.k_z = np.append(kz, kz)

            # Save this for future reference (we'll be back)
            self._air_anallos[(period)] = an
            return an
















def dec_float_str(dec_float):
    """ Convert float with decimal point into string with '_' in place of '.' """
    string = str(dec_float)
    fmt_string = string.replace('.','_')
    return fmt_string



def calculate_ff(square, d, a1, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0, a8=0, a9=0, a10=0,
    a11=0, a12=0, a13=0, a14=0, a15=0, a16=0, el1 = 0):

    if square == False:
        ff = np.pi*((a1/2)**2*np.sqrt(1-el1) + (a2/2)**2 + (a3/2)**2 + (a4/2)**2 + (a5/2)**2 + (a6/2)**2 + 
            (a7/2)**2 + (a8/2)**2 + (a9/2)**2 + (a10/2)**2 + (a11/2)**2 + (a12/2)**2 + (a13/2)**2 + 
            (a14/2)**2 + (a15/2)**2 + (a16/2)**2)/(d)**2
    else:
        ff = ((a1)**2 + (a2)**2 + (a3)**2 + (a4)**2 + (a5)**2 + (a6)**2 + (a7)**2 + (a8)**2 + (a9)**2
            + (a10)**2 + (a11)**2 + (a12)**2 + (a13)**2 + (a14)**2 + (a15)**2 + (a16)**2)/(d)**2 
    return ff