"""
    objects.py is a subroutine of EMUstack that contains the NanoStruct,
    ThinFilm and Light objects. These represent the properties of a
    structured layer, a homogeneous layer and the incident light
    respectively.

    Copyright (C) 2015  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

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
from __future__ import print_function

import os
import numpy as np
import paths
import materials
from mode_calcs import Simmo, Anallo
from fortran import EMUstack

msh_location = paths.msh_path
template_location = paths.template_path

# Acknowledgements
print('\n##################################################################\n'\
    + 'EMUstack is brought to you by Bjorn Sturmberg, Kokou Dossou, \n' \
    + 'Felix Lawrence & Lindsay Botton, with support from CUDOS & ARENA\n' \
    + 'Starting EMUstack calculation ...\n' + \
      '##################################################################\n')


class NanoStruct(object):
    """ Represents a structured layer.

        Args:
            periodicity  (str): Either 1D or 2D structure '1D_array', '2D_array'.

            period  (float): The period of the unit cell in nanometers.

            diameter1  (float): The diameter of the inclusion in nm.

        Keyword Args:
            period_y  (float): The period of the unit cell in the y-direction.\
                If None, period_y = period.

            inc_shape  (str): Shape of inclusions that have template mesh, \
                currently; 'circle', 'ellipse', 'square', 'ring', 'SRR',
                'dimer', 'square_dimer', 'strip_circle', 'strip_square',
                'rectangle', 'rectangle_shell', 'square_dimer_shell',
                'cross', 'cross_shell', 'L' .

            is_hex  (bool): Simulating a hexagonal lattice, using a rect unitcell?

            ellipticity  (float): If != 0, inclusion has given ellipticity, \
                with b = diameter, a = diameter-ellipticity * diameter. \
                NOTE: only implemented for a single inclusion.

            len_vertical  (float): Vertical length of split ring resonator \
                (if inc_shape = 'SRR').

            len_horizontal  (float): Horizontal length of split ring resonator\
                (if inc_shape = 'SRR').

            diameter2-16  (float): The diameters of further inclusions in nm. \
                Implemented up to diameter6 for 1D_arrays.

            gap (float): The dimer gap in nm. \
                (if inc_shape = 'dimer' or 'square_dimer').

            smooth (float): smoothness of square_dimer angles, between 0 (sharp). \
                and 1 (circle).
                (if inc_shape = 'square_dimer' or inc_shape = 'rectangle').

            t (float): shell thickness. \
                (if inc_shape = 'every _shell target').

            inclusion_a  : A :Material: instance for first inclusion, \
                specified as dispersive refractive index (eg. materials.Si_c) \
                or nondispersive complex number (eg. Material(1.0 + 0.0j)).

            inclusion_b  : A :Material: instance for the second \
                inclusion medium.

            inclusion_c  : A :Material: instance for the third \
                inclusion medium.

            inclusion_d  : A :Material: instance for the fourth \
                inclusion medium.

            inclusion_e  : A :Material: instance for the fifth \
                inclusion medium.

            background  : A :Material: instance for the background medium.

            loss  (bool): If False, Im(n) = 0, if True n as in \
                :Material: instance.

            height_nm  (float): The thickness of the layer in nm or \
                'semi_inf' for a semi-infinite layer.

            hyperbolic  (bool): If True FEM looks for Eigenvalues around \
                n**2 * k_0**2 rather than the regular \
                n**2 * k_0**2 - alpha**2 - beta**2.

            world_1d  (bool): Does the rest of the stack have exclusively 1D \
                periodic structures and homogeneous layers? \
                If True we use the set of 1D diffraction order PWs.\
                Defaults to True for '1D_array', and False for '2D_array'.

            ff  (float): The fill fraction of the inclusions. If non-zero, \
                the specified diameters are overwritten s.t. given ff is \
                achieved, otherwise ff is calculated from parameters and \
                stored in self.ff.

            ff_rand  (bool): If True, diameters overwritten with random \
                diameters, s.t. the ff is as assigned. Must provide non-zero \
                dummy diameters.

            posx  (float): Shift NWs laterally towards center (each other), \
                posx is a fraction of the distance possible before NWs touch.

            posy  (float): Shift NWs vertically towards center (each other), \
                posx is a fraction of the distance possible before NWs touch.

            small_space  (float): Only for 1D_arrays with 2 interleaved \
                inclusions. Sets distance between edges of inclusions. \
                By default (d_in_nm - diameter1 - diameter2) / 2. \
                The smaller distance is on the, which left of center \
                (inclusion_a remains centered).

            edge_spacing  (bool): For 1D_array with >= 3 inclusions. Space \
                inclusion surfaces by equal separations. Else their centers \
                will be equally spaced.

            split_touching_incs  (bool): For 1D_array with > 1 inclusions. \
                Arrange inclusions with touching edges, with the \
                aggregate centered in the unit cell.

            make_mesh_now  (bool): If True, program creates a FEM mesh with \
                provided :NanoStruct: parameters. If False, must provide \
                mesh_file name of existing .mail that will be run despite \
                :NanoStruct: parameters.

            force_mesh  (bool): If True, a new mesh is created despite \
                existence of mesh with same parameter. This is used to make \
                mesh with equal period etc. but different lc refinement.

            mesh_file  (str): If using a set premade mesh give its name \
                including .mail if 2D_array (eg. 600_60.mail), or .txt if \
                1D_array. It must be located in backend/fortran/msh/

            lc_bkg  (float): Length constant of meshing of background medium \
                (smaller = finer mesh)

            lc2  (float): factor by which lc_bkg will be reduced on inclusion \
                surfaces; lc_surface = cl_bkg / lc2.

            lc3-6'  (float): factor by which lc_bkg will be reduced at center \
                of inclusions.

            plotting_fields  (bool): Unless set to true field data deleted.\
                Also plots modes (ie. FEM solutions) in gmsh format. \
                Plots epsilon*|E|^2 & choice of real/imag/abs of \
                x,y,z components & field vectors. Fields are saved as gmsh \
                files, but can be converted by running the .geo file found in \
                Bloch_fields/PNG/

            plot_real  (bool): Choose to plot real part of modal fields.

            plot_imag  (bool): Choose to plot imaginary part of modal fields.

            plot_abs  (bool): Choose to plot absolute value of modal fields.

            plt_msh  (bool): Save a plot of the 1D array geometry.
    """

    def __init__(self,
                 periodicity,
                 period,
                 diameter1,
                 period_y=None,
                 inc_shape='circle',
                 is_hex=False,
                 ellipticity=0.0,
                 ff=0,
                 ff_rand=False,
                 small_space=None,
                 edge_spacing=False,
                 split_touching_incs=False,
                 len_vertical=0,
                 len_horizontal=0,
                 background=materials.Material(1.0 + 0.0j),
                 inclusion_a=materials.Material(1.0 + 0.0j),
                 inclusion_b=materials.Material(1.0 + 0.0j),
                 inclusion_c=materials.Material(1.0 + 0.0j),
                 inclusion_d=materials.Material(1.0 + 0.0j),
                 inclusion_e=materials.Material(1.0 + 0.0j),
                 inclusion_f=materials.Material(1.0 + 0.0j),
                 loss=True,
                 height_nm=100.0,
                 diameter2=0,
                 diameter3=0,
                 diameter4=0,
                 diameter5=0,
                 diameter6=0,
                 diameter7=0,
                 diameter8=0,
                 diameter9=0,
                 diameter10=0,
                 diameter11=0,
                 diameter12=0,
                 diameter13=0,
                 diameter14=0,
                 diameter15=0,
                 diameter16=0,
                 gap=0,
                 smooth=0,
                 t=0,
                 hyperbolic=False,
                 world_1d=None,
                 posx=0,
                 posy=0,
                 xshift=None,
                 make_mesh_now=True,
                 force_mesh=True,
                 mesh_file='NEED_FILE.mail',
                 lc_bkg=0.09,
                 lc2=1.0,
                 lc3=1.0,
                 lc4=1.0,
                 lc5=1.0,
                 lc6=1.0,
                 plotting_fields=False,
                 plot_real=1,
                 plot_imag=0,
                 plot_abs=0,
                 plot_field_conc=False,
                 plt_msh=True):
        self.periodicity = periodicity
        self.period = float(period)
        self.diameter1 = diameter1
        if period_y is None:
            self.period_y = float(period)
        else:
            self.period_y = float(period_y)
        self.inc_shape = inc_shape
        self.is_hex = is_hex
        self.height_nm = height_nm
        self.background = background
        self.inclusion_a = inclusion_a
        self.inclusion_b = inclusion_b
        self.inclusion_c = inclusion_c
        self.inclusion_d = inclusion_d
        self.inclusion_e = inclusion_e
        self.inclusion_f = inclusion_f
        self.loss = loss
        self.hyperbolic = hyperbolic
        self.diameter2 = diameter2
        self.diameter3 = diameter3
        self.diameter4 = diameter4
        self.diameter5 = diameter5
        self.diameter6 = diameter6
        self.diameter7 = diameter7
        self.diameter8 = diameter8
        self.diameter9 = diameter9
        self.diameter10 = diameter10
        self.diameter11 = diameter11
        self.diameter12 = diameter12
        self.diameter13 = diameter13
        self.diameter14 = diameter14
        self.diameter15 = diameter15
        self.diameter16 = diameter16
        self.gap = gap
        self.smooth = smooth
        self.t = t
        self.len_vertical = len_vertical
        self.len_horizontal = len_horizontal
        self.ellipticity = ellipticity
        if ellipticity > 1.0:
            raise ValueError("ellipticity must be less than 1.0")
        if diameter3 != 0:
            self.nb_typ_el = 4
        elif diameter2 != 0:
            self.nb_typ_el = 3
        else:
            self.nb_typ_el = 2
        if ff == 0:
            if periodicity == '2D_array':
                                self.ff = calculate_ff(
                                    inc_shape, period, self.period_y, diameter1, diameter2,
                                    diameter3, diameter4, diameter5, diameter6, diameter7,
                                    diameter8, diameter9, diameter10, diameter11, diameter12,
                                    diameter13, diameter14, diameter15, diameter16,
                                    ellipticity)
            elif periodicity == '1D_array':
                self.ff = (diameter1 + diameter2) / period
        else:
            self.ff = ff
            if diameter2 != 0:
                self.diameter2 = 2 * ((ff * (period)**2) / np.pi -
                                      ((diameter1 / 2)**2))**0.5
            else:
                self.diameter1 = 2 * np.sqrt((ff * (period)**2) / np.pi)
        self.ff_rand = ff_rand
        if world_1d is None:
            if periodicity == '1D_array':
                self.world_1d = True
            if periodicity == '2D_array':
                self.world_1d = False
        else:
            self.world_1d = world_1d
        self.posx = posx
        self.posy = posy
        self.lc = lc_bkg
        self.lc2 = lc2
        self.lc3 = lc3
        self.lc4 = lc4
        self.lc5 = lc5
        self.lc6 = lc6
        self.force_mesh = force_mesh
        self.small_space = small_space
        self.edge_spacing = edge_spacing
        self.split_touching_incs = split_touching_incs
        self.plt_msh = plt_msh
        if make_mesh_now is True:
            self.make_mesh()
        else:
            self.mesh_file = mesh_file
        if plotting_fields is True:
            self.plotting_fields = 1
            if periodicity == '2D_array':
                if not os.path.exists("Bloch_fields"):
                    os.mkdir("Bloch_fields")
                if not os.path.exists("Bloch_fields/PDF"):
                    os.mkdir("Bloch_fields/PDF")
        else:
            self.plotting_fields = 0
        self.plot_real = plot_real
        self.plot_imag = plot_imag
        self.plot_abs = plot_abs
        self.plot_field_conc = plot_field_conc
        self.xshift = xshift

    def make_mesh(self):
        if self.periodicity == '2D_array':
            if self.inc_shape in ['circle', 'ellipse', 'square']:
                if self.diameter10 > 0:
                    supercell = 16
                    msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s_%(diassss)s' % {
                   'd': dec_float_str(self.period),
                   'dy': dec_float_str(self.period_y),
                   'dia': dec_float_str(self.diameter1),
                   'dias': dec_float_str(self.diameter2),
                   'dias': dec_float_str(self.diameter2),
                   'diass': dec_float_str(self.diameter3),
                   'diasss': dec_float_str(self.diameter4),
                   'diassss': dec_float_str(self.diameter5)}
                elif self.diameter5 > 0:
                    supercell = 9
                    msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s_%(diassss)s' % {
                   'd': dec_float_str(self.period),
                   'dy': dec_float_str(self.period_y),
                   'dia': dec_float_str(self.diameter1),
                   'dias': dec_float_str(self.diameter2),
                   'diass': dec_float_str(self.diameter3),
                   'diasss': dec_float_str(self.diameter4),
                   'diassss': dec_float_str(self.diameter5)}
                elif self.diameter4 > 0:
                    supercell = 4
                    msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s' % {
                   'd': dec_float_str(self.period),
                   'dy': dec_float_str(self.period_y),
                   'dia': dec_float_str(self.diameter1),
                   'dias': dec_float_str(self.diameter2),
                   'diass': dec_float_str(self.diameter3),
                   'diasss': dec_float_str(self.diameter4)}
                elif self.diameter3 > 0:
                    supercell = 3
                    msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s' % {
                   'd': dec_float_str(self.period),
                   'dy': dec_float_str(self.period_y),
                   'dia': dec_float_str(self.diameter1),
                   'dias': dec_float_str(self.diameter2),
                   'diass': dec_float_str(self.diameter3)}
                elif self.diameter2 > 0:
                    supercell = 2
                    msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s' % {
                    'd': dec_float_str(self.period),
                   'dy': dec_float_str(self.period_y),
                   'dia': dec_float_str(self.diameter1),
                   'diameters': dec_float_str(self.diameter2)}
                elif self.diameter1 > 0:
                    supercell = 1
                    if self.is_hex is False:
                        msh_name = '%(d)s_%(dy)s_%(dia)s' % {
                                   'd': dec_float_str(self.period),
                                   'dy': dec_float_str(self.period_y),
                                   'dia': dec_float_str(self.diameter1)}
                    elif self.is_hex is True:
                        msh_name = 'hex_%(d)s_%(dy)s_%(dia)s' % {
                                   'd': dec_float_str(self.period),
                                   'dy': dec_float_str(self.period_y),
                                   'dia': dec_float_str(self.diameter1)}
                else:
                    raise ValueError("must have at least one cylinder of nonzero diameter.")

                if self.ellipticity != 0:
                    msh_name = msh_name + '_e_%(e)s' % {'e': dec_float_str(self.ellipticity),}
                if self.inc_shape == 'square':
                    msh_name = msh_name + '_sq'
                if self.posx != 0:
                    msh_name = msh_name + 'x%(e)s' % {'e': dec_float_str(self.posx),}
                if self.posy != 0:
                    msh_name = msh_name + 'y%(e)s' % {'e': dec_float_str(self.posy),}

                # for blah in range(1,101,1):
                #     print blah
                #     msh_name = 'random_u_%i' % blah
                #     self.mesh_file = msh_name + '.mail'
                # msh_name = 'design-last_17'
                if self.ff_rand is True:
                    import random
                    ff_tol = 0.0001
                    min_a = 50
                    max_a = (self.period/1.05)/np.sqrt(supercell)
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

                        test_ff = calculate_ff(self.inc_shape, self.period,
                                               self.period_y, rad_array[0],
                                               rad_array[1], rad_array[2],
                                               rad_array[3], rad_array[4],
                                               rad_array[5], rad_array[6],
                                               rad_array[7], rad_array[8],
                                               rad_array[9], rad_array[10],
                                               rad_array[11], rad_array[12],
                                               rad_array[13], rad_array[14],
                                               rad_array[15])
                        print(test_ff)
                        if supercell > 3:
                            self.diameter1 = rad_array[0]
                            self.diameter2 = rad_array[1]
                            self.diameter3 = rad_array[2]
                            self.diameter4 = rad_array[3]
                        if supercell > 4:
                            self.diameter5 = rad_array[4]
                            self.diameter6 = rad_array[5]
                            self.diameter7 = rad_array[6]
                            self.diameter8 = rad_array[7]
                            self.diameter9 = rad_array[8]
                        if supercell > 9:
                            self.diameter10 = rad_array[9]
                            self.diameter11 = rad_array[10]
                            self.diameter12 = rad_array[11]
                            self.diameter13 = rad_array[12]
                            self.diameter14 = rad_array[13]
                            self.diameter15 = rad_array[14]
                            self.diameter16 = rad_array[15]
                        test_ff = calculate_ff(self.inc_shape, self.period,
                                               self.period_y, rad_array[0],
                                               rad_array[1], rad_array[2],
                                               rad_array[3], rad_array[4],
                                               rad_array[5], rad_array[6],
                                               rad_array[7], rad_array[8],
                                               rad_array[9], rad_array[10],
                                               rad_array[11], rad_array[12],
                                               rad_array[13], rad_array[14],
                                               rad_array[15])


                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    if self.is_hex is False:
                        geo_tmp = open(template_location + '%s_msh_template.geo' % supercell, "r").read()
                    else:
                        geo_tmp = open(template_location + 'hex_msh_template.geo', "r").read()

                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;', "d_in_nm = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('ellipticity = 0;', "ellipticity = %f;" % self.ellipticity)
                    if self.inc_shape == 'square': geo = geo.replace('square = 0;', "square = 1;")
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


            elif self.inc_shape == 'SRR':
                msh_name = 'SRR_%(d)s_%(dy)s_%(lvert)s_%(lhori)s_%(dia)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'lvert': dec_float_str(self.len_vertical),
                           'lhori': dec_float_str(self.len_horizontal),
                           'dia': dec_float_str(self.diameter1)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + 'SRR_msh_template.geo', "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm  = 0;', "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('lvert_nm = 0;', "lvert_nm = %f;" % self.len_vertical)
                    geo = geo.replace('lhori_nm = 0;', "lhori_nm = %f;" % self.len_horizontal)
                    geo = geo.replace('width_nm = 0;', "width_nm = %f;" % self.diameter1)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)


            elif self.inc_shape == 'ring':
                msh_name = 'ring_%(d)s_%(dy)s_%(dia_out)s_%(dia_in)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'dia_out': dec_float_str(self.diameter1),
                           'dia_in': dec_float_str(self.diameter2)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + 'ring1_msh_template.geo', "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;', "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                    geo = geo.replace('xshift_nm = 0;', "xshift_nm = %f;" % self.xshift)

            elif self.inc_shape == 'dimer':
                msh_name = 'dimer_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(gap)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'd_one': dec_float_str(self.diameter1),
                           'd_two': dec_float_str(self.diameter2),
                           'gap': dec_float_str(self.gap)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + 'dimer1_msh_template.geo', "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;', "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('gap = 0;', "gap = %f;" % self.gap)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)

            elif self.inc_shape == 'square_dimer':
                msh_name = 'square_dimer_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(gap)s_%(smooth)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'd_one': dec_float_str(self.diameter1),
                           'd_two': dec_float_str(self.diameter2),
                           'gap': dec_float_str(self.gap),
                           'smooth': dec_float_str(self.smooth)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + 'square_dimer1_msh_template.geo', "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;', "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('gap = 0;', "gap = %f;" % self.gap)
                    geo = geo.replace('smooth = 0;', "smooth = %f;" % self.smooth)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)

            elif self.inc_shape == 'square_shell_dimer':
                msh_name = 'square_shell_dimer_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(d_three)s_%(d_four)s_%(gap)s_%(smooth)s_%(t)s' % {
                    'd': dec_float_str(self.period),
                    'dy': dec_float_str(self.period_y),
                    'd_one': dec_float_str(self.diameter1),
                    'd_two': dec_float_str(self.diameter2),
                    'd_three': dec_float_str(self.diameter3),
                    'd_four': dec_float_str(self.diameter4),
                    'gap': dec_float_str(self.gap),
                    'smooth': dec_float_str(self.smooth),
                    't': dec_float_str(self.t)
                }
                if not os.path.exists(msh_location + msh_name +
                                      '.mail') or self.force_mesh is True:
                    geo_tmp = open(
                        template_location + 'square_shell_dimer1_msh_template.geo',
                        "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;',
                                      "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;',
                                      "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('b1 = 0;', "b1 = %f;" % self.diameter3)
                    geo = geo.replace('b2 = 0;', "b2 = %f;" % self.diameter4)
                    geo = geo.replace('gap = 0;', "gap = %f;" % self.gap)
                    geo = geo.replace('smooth = 0;',
                                      "smooth = %f;" % self.smooth)
                    geo = geo.replace('t = 0;', "t = %f;" % self.t)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
            elif self.inc_shape == 'rectangle':
                msh_name = 'rect_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(smooth)s' % {
                    'd': dec_float_str(self.period),
                    'dy': dec_float_str(self.period_y),
                    'd_one': dec_float_str(self.diameter1),
                    'd_two': dec_float_str(self.diameter2),
                    'smooth': dec_float_str(self.smooth)
                }
                if not os.path.exists(msh_location + msh_name +
                                      '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + 'rect1_msh_template.geo',
                                   "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;',
                                      "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;',
                                      "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('smooth = 0;',
                                      "smooth = %f;" % self.smooth)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
            elif self.inc_shape == 'rectangle_shell':
                msh_name = 'rect_shell_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(smooth)s_%(t)s' % {
                    'd': dec_float_str(self.period),
                    'dy': dec_float_str(self.period_y),
                    'd_one': dec_float_str(self.diameter1),
                    'd_two': dec_float_str(self.diameter2),
                    'smooth': dec_float_str(self.smooth),
                    't': dec_float_str(self.t)
                }
                if not os.path.exists(msh_location + msh_name +
                                      '.mail') or self.force_mesh is True:
                    geo_tmp = open(
                        template_location + 'rect_shell1_msh_template.geo',
                        "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;',
                                      "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;',
                                      "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('b1 = 0;', "b1 = %f;" % self.diameter2)
                    geo = geo.replace('smooth = 0;',
                                      "smooth = %f;" % self.smooth)
                    geo = geo.replace('t = 0;', "t = %f;" % self.t)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
            elif self.inc_shape == 'cross':
                msh_name = 'cross_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(smooth)s' % {
                    'd': dec_float_str(self.period),
                    'dy': dec_float_str(self.period_y),
                    'd_one': dec_float_str(self.diameter1),
                    'd_two': dec_float_str(self.diameter2),
                    'smooth': dec_float_str(self.smooth)
                }
                if not os.path.exists(msh_location + msh_name +
                                      '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + 'cross1_msh_template.geo',
                                   "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;',
                                      "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;',
                                      "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('smooth = 0;',
                                      "smooth = %f;" % self.smooth)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
            elif self.inc_shape == 'cross_shell':
                msh_name = 'cross_shell_%(d)s_%(dy)s_%(d_one)s_%(d_two)s_%(smooth)s_%(t)s' % {
                    'd': dec_float_str(self.period),
                    'dy': dec_float_str(self.period_y),
                    'd_one': dec_float_str(self.diameter1),
                    'd_two': dec_float_str(self.diameter2),
                    'smooth': dec_float_str(self.smooth),
                    't': dec_float_str(self.t)
                }
                if not os.path.exists(msh_location + msh_name +
                                      '.mail') or self.force_mesh is True:
                    geo_tmp = open(
                        template_location + 'cross_shell1_msh_template.geo',
                        "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;',
                                      "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;',
                                      "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('smooth = 0;',
                                      "smooth = %f;" % self.smooth)
                    geo = geo.replace('t = 0;', "t = %f;" % self.t)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)

            elif self.inc_shape == 'L':
                msh_name = 'L_%(d)s_%(dy)s_%(L)s_%(W)s_%(r)s' % {
                    'd': dec_float_str(self.period),
                    'dy': dec_float_str(self.period_y),
                    'L': dec_float_str(self.diameter1),
                    'W': dec_float_str(self.diameter2),
                    'r': dec_float_str(self.smooth),

                }
                if not os.path.exists(msh_location + msh_name +
                                      '.mail') or self.force_mesh is True:
                    geo_tmp = open(
                        template + 'L_msh_template.geo',
                        "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;',
                                      "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;',
                                      "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('L_nm = 0;', "L_nm = %f;" % self.diameter1)
                    geo = geo.replace('W_nm = 0;', "W_nm = %f;" % self.diameter2)
                    geo = geo.replace('r = 0;',
                                      "r = %f;" % self.smooth)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
            elif self.inc_shape == 'strip_circle':
                msh_name = 'strip_circle_%(d)s_%(dy)s_%(d_one)s_%(d_two)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'd_one': dec_float_str(self.diameter1),
                           'd_two': dec_float_str(self.diameter2)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + '1_strip_msh_template.geo', "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;', "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('strip = 0;', "strip = %f;" % self.diameter2)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)

            elif self.inc_shape == 'strip_square':
                msh_name = 'strip_square_%(d)s_%(dy)s_%(d_one)s_%(d_two)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'd_one': dec_float_str(self.diameter1),
                           'd_two': dec_float_str(self.diameter2)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(template_location + '1_strip_msh_template.geo', "r").read()
                    geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                    geo = geo.replace('d_in_nm = 0;', "d_in_nm  = %f;" % self.period)
                    geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                    geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                    geo = geo.replace('strip = 0;', "strip = %f;" % self.diameter2)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                    geo = geo.replace('square = 0;', "square = 1;")

            else:
                raise NotImplementedError("\n Selected inc_shape = '%s' \n \
                is not currently implemented. Please make a mesh with gmsh, & \n \
                consider contributing this to EMUstack via github." % self.inc_shape)

            self.mesh_file = msh_name + '.mail'
            if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                open(msh_location + msh_name + '.geo', "w").write(geo)
                EMUstack.conv_gmsh(msh_location+msh_name)

            # # Automatically show created mesh in gmsh.
            # gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.msh'
            # os.system(gmsh_cmd)
            # gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.geo'
            # os.system(gmsh_cmd)


        elif self.periodicity == '1D_array':
            # Unit cell length normalized to unity
            x_min = 0.0
            x_max = 1.0
            # Mesh elements and points
            nel = int(np.round(1.0/self.lc))
            npt = 2 * nel + 1
            delta_x = (x_max - x_min) / nel
            # Coordinate and type of the nodes
            el_list = range(1,nel+1)
            table_nod = np.zeros((3,nel+1))
            type_el = np.zeros(nel+1)
            ls_x = np.zeros(npt+1)

            for i_el in el_list:
                x = x_min + (i_el-1) * delta_x
                ls_x[2*i_el-1] = x
                ls_x[2*i_el] = x + delta_x / 2.0
            # End-points
            x = x_min + i_el * delta_x
            ls_x[2*i_el+1] = x
            # Connectivity table
            for i_el in el_list:
                table_nod[0, i_el] = 2*i_el-1
                table_nod[1, i_el] = 2*i_el+1
                table_nod[2, i_el] = 2*i_el  # Mid-node

            if self.diameter6 > 0:
                msh_name = '%(d)s_%(di)s_%(dis)s_%(diss)s_%(disss)s_%(dissss)s_%(disssss)s' % {
                'd': dec_float_str(self.period), 'di': dec_float_str(self.diameter1),
                'dis': dec_float_str(self.diameter2), 'diss': dec_float_str(self.diameter3),
                'disss': dec_float_str(self.diameter4), 'dissss': dec_float_str(self.diameter5),
                'disssss': dec_float_str(self.diameter6)}
                # End-points of the elements
                rad_1 = self.diameter1/(2.0*self.period)
                rad_2 = self.diameter2/(2.0*self.period)
                rad_3 = self.diameter3/(2.0*self.period)
                rad_4 = self.diameter4/(2.0*self.period)
                rad_5 = self.diameter5/(2.0*self.period)
                rad_6 = self.diameter6/(2.0*self.period)
                if self.edge_spacing is True:
                    i_d = 2.0*(0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 - rad_6)/6.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif x_1 <= 0.5 - i_d - rad_1 and 0.5 - i_d - rad_1 - 2.0*rad_2 <= x_1 \
                        and  x_2 <= 0.5 - i_d - rad_1 and 0.5 - i_d - rad_1 - 2.0*rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif x_1 <= 0.5 + i_d + 2.0*rad_3 + rad_1 and 0.5 + i_d + rad_1 <= x_1 \
                        and  x_2 <= 0.5 + i_d + 2.0*rad_3 + rad_1 and 0.5 + i_d + rad_1 <= x_2:
                            type_el[i_el] = 4
                        elif x_1 <= 0.5 - 2.0*i_d - 2.0*rad_2 - rad_1 and 0.5 - 2.0*i_d - 2.0*rad_2 - 2.0*rad_4 - rad_1 <= x_1 \
                        and  x_2 <= 0.5 - 2.0*i_d - 2.0*rad_2 - rad_1 and 0.5 - 2.0*i_d - 2.0*rad_2 - 2.0*rad_4 - rad_1 <= x_2:
                            type_el[i_el] = 5
                        elif 0.5 + 2.0*i_d + 2.0*rad_3 + rad_1 <= x_1 and x_1 <= 0.5 + 2.0*i_d + 2.0*rad_3 + 2.0*rad_5 + rad_1 \
                        and  0.5 + 2.0*i_d + 2.0*rad_3 + rad_1 <= x_2 and x_2 <= 0.5 + 2.0*i_d + 2.0*rad_3 + 2.0*rad_5 + rad_1:
                            type_el[i_el] = 6
                        elif x_1 <= 0.5 - 3.0*i_d - rad_1 - 2.0*rad_2 - 2.0*rad_4 and x_1 >= 0.5 - 3.0*i_d - rad_1 - 2.0*rad_2 - 2.0*rad_4 - rad_5\
                        and  x_2 <= 0.5 - 3.0*i_d - rad_1 - 2.0*rad_2 - 2.0*rad_4 and x_2 >= 0.5 - 3.0*i_d - rad_1 - 2.0*rad_2 - 2.0*rad_4 - rad_5:
                            type_el[i_el] = 6
                        elif x_1 >= 0.5 + 3.0*i_d + rad_1 + 2.0*rad_3 + 2.0*rad_5 and x_1 <= 0.5 + 3.0*i_d + rad_1 + 2.0*rad_3 + 2.0*rad_5 + 2.0*rad_6\
                        and  x_2 >= 0.5 + 3.0*i_d + rad_1 + 2.0*rad_3 + 2.0*rad_5 and x_2 <= 0.5 + 3.0*i_d + rad_1 + 2.0*rad_3 + 2.0*rad_5 + 2.0*rad_6:
                            type_el[i_el] = 7
                        else:
                            type_el[i_el] = 1
                elif self.split_touching_incs is True:
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif  x_1 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_1 \
                        and x_2 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 and 0.5 + rad_1 + rad_2 + rad_3 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 and 0.5 + rad_1 + rad_2 + rad_3 <= x_2:
                            type_el[i_el] = 5
                        elif  x_1 <= 0.5 - rad_1 - rad_2 - rad_3 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 - rad_3 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 <= x_2:
                            type_el[i_el] = 5
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 + rad_5 and 0.5 + rad_1 + rad_2 + rad_3 + rad_4 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 + rad_5 and 0.5 + rad_1 + rad_2 + rad_3 + rad_4 <= x_2:
                            type_el[i_el] = 6
                        elif  x_1 <= 0.5 - rad_1 - rad_2 - rad_3 - rad_4 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 - rad_3 - rad_4 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 <= x_2:
                            type_el[i_el] = 6
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 + rad_5 + rad_6 and 0.5 + rad_1 + rad_2 + rad_3 + rad_4 +rad_5 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 + rad_5 + rad_6 and 0.5 + rad_1 + rad_2 + rad_3 + rad_4 +rad_5 <= x_2:
                            type_el[i_el] = 7
                        elif  x_1 <= 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 -rad_6 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 - rad_3 - rad_4 -rad_5 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 -rad_6 <= x_2:
                            type_el[i_el] = 7
                        else:
                            type_el[i_el] = 1
                else:
                    i_d = 1.0/6.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif 0.5 - i_d - rad_2 <= x_1 and x_1 <= 0.5 - i_d + rad_2 \
                        and  0.5 - i_d - rad_2 <= x_2 and x_2 <= 0.5 - i_d + rad_2:
                            type_el[i_el] = 3
                        elif 0.5 + i_d - rad_3 <= x_1 and x_1 <= 0.5 + i_d + rad_3 \
                        and  0.5 + i_d - rad_3 <= x_2 and x_2 <= 0.5 + i_d + rad_3:
                            type_el[i_el] = 4
                        elif 0.5 - 2.0*i_d - rad_4 <= x_1 and x_1 <= 0.5 - 2.0*i_d + rad_4 \
                        and  0.5 - 2.0*i_d - rad_4 <= x_2 and x_2 <= 0.5 - 2.0*i_d + rad_4:
                            type_el[i_el] = 5
                        elif 0.5 + 2.0*i_d - rad_5 <= x_1 and x_1 <= 0.5 + 2.0*i_d + rad_5 \
                        and  0.5 + 2.0*i_d - rad_5 <= x_2 and x_2 <= 0.5 + 2.0*i_d + rad_5:
                            type_el[i_el] = 6
                        elif x_1 <= rad_6 and x_2 <= rad_6:
                            type_el[i_el] = 7
                        elif x_1 >= 1.0 - rad_6 and x_2 >= 1.0 - rad_6:
                            type_el[i_el] = 7
                        else:
                            type_el[i_el] = 1
            elif self.diameter5 > 0:
                msh_name = '%(d)s_%(di)s_%(dis)s_%(diss)s_%(disss)s_%(dissss)s' % {
                           'd': dec_float_str(self.period),
                           'di': dec_float_str(self.diameter1),
                           'dis': dec_float_str(self.diameter2),
                           'diss': dec_float_str(self.diameter3),
                           'disss': dec_float_str(self.diameter4),
                           'dissss': dec_float_str(self.diameter5)}
                # End-points of the elements
                rad_1 = self.diameter1/(2.0*self.period)
                rad_2 = self.diameter2/(2.0*self.period)
                rad_3 = self.diameter3/(2.0*self.period)
                rad_4 = self.diameter4/(2.0*self.period)
                rad_5 = self.diameter5/(2.0*self.period)
                if self.edge_spacing is True:
                    i_d = 2.0*(0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5)/5.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif x_1 <= 0.5 - i_d - rad_1 and 0.5 - i_d - rad_1 - 2.0*rad_2 <= x_1 \
                        and x_2 <= 0.5 - i_d - rad_1 and 0.5 - i_d - rad_1 - 2.0*rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif x_1 <= 0.5 + i_d + 2.0*rad_3 + rad_1 and 0.5 + i_d + rad_1 <= x_1 \
                        and x_2 <= 0.5 + i_d + 2.0*rad_3 + rad_1 and 0.5 + i_d + rad_1 <= x_2:
                            type_el[i_el] = 4
                        elif x_1 <= 0.5 - 2.0*i_d - 2.0*rad_2 - rad_1 and 0.5 - 2.0*i_d - 2.0*rad_2 - 2.0*rad_4 - rad_1 <= x_1 \
                        and x_2 <= 0.5 - 2.0*i_d - 2.0*rad_2 - rad_1 and 0.5 - 2.0*i_d - 2.0*rad_2 - 2.0*rad_4 - rad_1 <= x_2:
                            type_el[i_el] = 5
                        elif 0.5 + 2.0*i_d + 2.0*rad_3 + rad_1 <= x_1 and x_1 <= 0.5 + 2.0*i_d + 2.0*rad_3 + 2.0*rad_5 + rad_1 \
                        and 0.5 + 2.0*i_d + 2.0*rad_3 + rad_1 <= x_2 and x_2 <= 0.5 + 2.0*i_d + 2.0*rad_3 + 2.0*rad_5 + rad_1:
                            type_el[i_el] = 4
                        else:
                            type_el[i_el] = 1
                elif self.split_touching_incs is True:
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif  x_1 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_1 \
                        and x_2 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 and 0.5 + rad_1 + rad_2 + rad_3 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 and 0.5 + rad_1 + rad_2 + rad_3 <= x_2:
                            type_el[i_el] = 5
                        elif  x_1 <= 0.5 - rad_1 - rad_2 - rad_3 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 - rad_3 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 <= x_2:
                            type_el[i_el] = 5
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 + rad_5 and 0.5 + rad_1 + rad_2 + rad_3 + rad_4 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 + rad_5 and 0.5 + rad_1 + rad_2 + rad_3 + rad_4 <= x_2:
                            type_el[i_el] = 6
                        elif  x_1 <= 0.5 - rad_1 - rad_2 - rad_3 - rad_4 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 - rad_3 - rad_4 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 - rad_5 <= x_2:
                            type_el[i_el] = 6
                        else:
                            type_el[i_el] = 1
                else:
                    i_d = 1.0/5.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif 0.5 - i_d - rad_2 <= x_1 and x_1 <= 0.5 - i_d + rad_2 \
                        and 0.5 - i_d - rad_2 <= x_2 and x_2 <= 0.5 - i_d + rad_2:
                            type_el[i_el] = 3
                        elif 0.5 + i_d - rad_3 <= x_1 and x_1 <= 0.5 + i_d + rad_3 \
                        and 0.5 + i_d - rad_3 <= x_2 and x_2 <= 0.5 + i_d + rad_3:
                            type_el[i_el] = 4
                        elif 0.5 - 2.0*i_d - rad_4 <= x_1 and x_1 <= 0.5 - 2.0*i_d + rad_4 \
                        and 0.5 - 2.0*i_d - rad_4 <= x_2 and x_2 <= 0.5 - 2.0*i_d + rad_4:
                            type_el[i_el] = 5
                        elif 0.5 + 2.0*i_d - rad_5 <= x_1 and x_1 <= 0.5 + 2.0*i_d + rad_5 \
                        and 0.5 + 2.0*i_d - rad_5 <= x_2 and x_2 <= 0.5 + 2.0*i_d + rad_5:
                            type_el[i_el] = 6
                        else:
                            type_el[i_el] = 1
            elif self.diameter4 > 0:
                msh_name = '%(d)s_%(di)s_%(dis)s_%(diss)s_%(disss)s' % {
                           'd': dec_float_str(self.period),
                           'di': dec_float_str(self.diameter1),
                           'dis': dec_float_str(self.diameter2),
                           'diss': dec_float_str(self.diameter3),
                           'disss': dec_float_str(self.diameter4)}
                # End-points of the elements
                rad_1 = self.diameter1/(2.0*self.period)
                rad_2 = self.diameter2/(2.0*self.period)
                rad_3 = self.diameter3/(2.0*self.period)
                rad_4 = self.diameter4/(2.0*self.period)
                if self.edge_spacing is True:
                    i_d = 2.0*(0.5 - rad_1 - rad_2 - rad_3 - rad_4)/4.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif x_1 <= 0.5 - i_d - rad_1 and 0.5 - i_d - rad_1 - 2.0*rad_2 <= x_1 \
                        and x_2 <= 0.5 - i_d - rad_1 and 0.5 - i_d - rad_1 - 2.0*rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif 0.5 + i_d + rad_1 <= x_1 and x_1 <= 0.5 + i_d + rad_1 + 2.0*rad_3 \
                        and 0.5 + i_d + rad_1 <= x_2 and x_2 <= 0.5 + i_d + rad_1 + 2.0*rad_3:
                            type_el[i_el] = 4
                        elif x_1 >= 0.5 + 2.0*i_d + rad_1 + 2.0*rad_3 \
                        and x_2 >= 0.5 + 2.0*i_d + rad_1 + 2.0*rad_3:
                            type_el[i_el] = 4
                        elif x_1 <= 0.5 - 2.0*i_d - rad_1 - 2.0*rad_2 \
                        and x_2 <= 0.5 - 2.0*i_d - rad_1 - 2.0*rad_2:
                            type_el[i_el] = 3
                        else:
                            type_el[i_el] = 1
                elif self.split_touching_incs is True:
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif  x_1 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_1 \
                        and x_2 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 and 0.5 + rad_1 + rad_2 + rad_3 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 + rad_4 and 0.5 + rad_1 + rad_2 + rad_3 <= x_2:
                            type_el[i_el] = 5
                        elif  x_1 <= 0.5 - rad_1 - rad_2 - rad_3 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 - rad_3 and 0.5 - rad_1 - rad_2 - rad_3 - rad_4 <= x_2:
                            type_el[i_el] = 5
                        else:
                            type_el[i_el] = 1
                else:
                    i_d = 1.0/4.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif 0.5 - i_d - rad_2 <= x_1 and x_1 <= 0.5 - i_d + rad_2 \
                        and  0.5 - i_d - rad_2 <= x_2 and x_2 <= 0.5 - i_d + rad_2:
                            type_el[i_el] = 3
                        elif 0.5 + i_d - rad_3 <= x_1 and x_1 <= 0.5 + i_d + rad_3 \
                        and  0.5 + i_d - rad_3 <= x_2 and x_2 <= 0.5 + i_d + rad_3:
                            type_el[i_el] = 4
                        elif x_1 <= rad_4 and x_2 <= rad_4:
                            type_el[i_el] = 5
                        elif x_1 >= 1.0 - rad_4 and x_2 >= 1.0 - rad_4:
                            type_el[i_el] = 5
                        else:
                            type_el[i_el] = 1
            elif self.diameter3 > 0:
                msh_name = '%(d)s_%(di)s_%(dis)s_%(diss)s' % {
                           'd': dec_float_str(self.period),
                           'di': dec_float_str(self.diameter1),
                           'dis': dec_float_str(self.diameter2),
                           'diss': dec_float_str(self.diameter3)}
                # End-points of the elements
                rad_1 = self.diameter1/(2.0*self.period)
                rad_2 = self.diameter2/(2.0*self.period)
                rad_3 = self.diameter3/(2.0*self.period)
                if self.edge_spacing is True:
                    i_d = (1.0 - self.diameter1 - self.diameter2 - self.diameter3)/3.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        # inclusion 1
                        if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        # inclusion 2
                        elif 0.5 - i_d - 2.0*rad_2 - rad_1 <= x_1 and x_1 <= 0.5 - i_d - rad_1 \
                        and  0.5 - i_d - 2.0*rad_2 - rad_1 <= x_2 and x_2 <= 0.5 - i_d - rad_1:
                            type_el[i_el] = 3
                        # inclusion 3
                        elif x_1 <= 0.5 + i_d + 2.0*rad_3 + rad_1 and 0.5 + i_d + rad_1 <= x_1 \
                        and  x_2 <= 0.5 + i_d + 2.0*rad_3 + rad_1 and 0.5 + i_d + rad_1 <= x_2:
                            type_el[i_el] = 4
                        else:
                            type_el[i_el] = 1
                elif self.split_touching_incs is True:
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif  x_1 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_1 \
                        and x_2 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 + rad_3 and 0.5 + rad_1 + rad_2 <= x_2:
                            type_el[i_el] = 4
                        elif  x_1 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_1 \
                        and x_2 <= 0.5 - rad_1 - rad_2 and 0.5 - rad_1 - rad_2 - rad_3 <= x_2:
                            type_el[i_el] = 4
                        else:
                            type_el[i_el] = 1
                else:
                    i_d = 1.0/3.0
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif 0.5 - i_d - rad_2 <= x_1 and x_1 <= 0.5 - i_d + rad_2 \
                        and  0.5 - i_d - rad_2 <= x_2 and x_2 <= 0.5 - i_d + rad_2:
                            type_el[i_el] = 3
                        elif 0.5 + i_d - rad_3 <= x_1 and x_1 <= 0.5 + i_d + rad_3 \
                        and  0.5 + i_d - rad_3 <= x_2 and x_2 <= 0.5 + i_d + rad_3:
                            type_el[i_el] = 4
                        else:
                            type_el[i_el] = 1
            elif self.diameter2 > 0:
                msh_name = '1D_%(d)s_%(diameter)s_%(diameters)s' % {
                           'd': dec_float_str(self.period),
                           'diameter': dec_float_str(self.diameter1),
                           'diameters': dec_float_str(self.diameter2)}
                # End-points of the elements
                rad_1 = self.diameter1/(2.0*self.period)
                rad_2 = self.diameter2/(2.0*self.period)
                if self.split_touching_incs is True:
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if  x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif  x_1 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 + rad_2 and 0.5 + rad_1 <= x_2:
                            type_el[i_el] = 3
                        elif  x_1 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_1 \
                        and x_2 <= 0.5 - rad_1 and 0.5 - rad_1 - rad_2 <= x_2:
                            type_el[i_el] = 3
                        else:
                            type_el[i_el] = 1
                else:
                    if self.small_space is None:
                        small_space = large_d = 0.5 - rad_1 - rad_2
                    else:
                        small_space = self.small_space
                        large_d = 1.0 - small_space - (2*rad_1) - (2*rad_2)
                    for i_el in el_list:
                        x_1 = ls_x[2*i_el-1]
                        x_2 = ls_x[2*i_el+1]
                        if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                        and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                            type_el[i_el] = 2
                        elif x_1 <= 0.5-rad_1-small_space and x_2 <= 0.5-rad_1-small_space:
                            type_el[i_el] = 3
                        elif x_1 >= 0.5+large_d+rad_1 and x_2 >= 0.5+large_d+rad_1:
                            type_el[i_el] = 3
                        else:
                            type_el[i_el] = 1
            elif self.diameter1 > 0:
                msh_name  =  '1D_%(d)s_%(diameter)s' % {'d' : dec_float_str(self.period),
                    'diameter' : dec_float_str(self.diameter1)}
                # End-points of the elements
                rad_1 = self.diameter1/(2.0*self.period)
                for i_el in el_list:
                    x_1 = ls_x[2*i_el-1]
                    x_2 = ls_x[2*i_el+1]
                    if x_1 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_1 \
                    and x_2 <= 0.5 + rad_1 and 0.5 - rad_1 <= x_2:
                        type_el[i_el] = 2
                    else:
                        type_el[i_el] = 1
            else:
                raise ValueError("Must have at least one grating of nonzero width.")

            # Store useful quantities as property of the object.
            self.n_msh_el = nel
            self.n_msh_pts = npt
            self.table_nod = table_nod[:,1:]
            # self.type_el = type_el[1:]
            self.type_el = type_el[1:]
            self.x_arr = ls_x[1:]
            self.mesh_file = msh_name

            if self.plt_msh is True:
                import matplotlib
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax1 = fig.add_subplot(1,1,1)
                ax1.plot(el_list,self.type_el)
                ax1.fill_between(el_list,self.type_el,0)
                ax1.set_xlim(el_list[0],el_list[-1])
                ax1.set_ylim(1,7)
                ax1.set_yticks([1,2,3,4,5,6,7])
                ax1.set_yticklabels(['bkg', 'inc_a', 'inc_b',
                    'inc_c', 'inc_d','inc_e', 'inc_f'])
                ax1.set_xlabel('Element Number')
                ax1.set_ylabel('Material Type')
                plt.savefig(msh_name, bbox_inches='tight')

            # Then clean up local variables.
            del nel, npt, table_nod, ls_x, type_el, el_list

        # Latency of old 1D grating meshed in 2D.

        # elif self.periodicity == '1D_array':
        #     if self.diameter2 > 0:
        #         supercell = 2
        #         msh_name  =  '1D_%(d)s_%(diameter)s_%(diameters)s' % {
        #        'd' : dec_float_str(self.period), 'diameter' : dec_float_str(self.diameter1),
        #        'diameters' : dec_float_str(self.diameter2)}
        #     elif self.diameter1 > 0:
        #         supercell = 1
        #         msh_name  =  '1D_%(d)s_%(diameter)s' % {'d' : dec_float_str(self.period),
        #             'diameter' : dec_float_str(self.diameter1)}
        #     else:
        #         raise ValueError, "must have at least one grating of nonzero width."

        #     self.mesh_file = msh_name + '.mail'


        #     if not os.path.exists(msh_location + msh_name + '.mail') or force_mesh == True:
        #         geo_tmp = open(msh_location + '1D_%s_msh_template.geo' % supercell, "r").read()
        #         geo = geo_tmp.replace('d_in_nm = 0;', "d_in_nm = %f;" % self.period)
        #         geo = geo.replace('w1 = 0;', "w1 = %f;" % self.diameter1)
        #         geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
        #         geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
        #         if supercell > 1:
        #             geo = geo.replace('w2 = 0;', "w2 = %f;" % self.diameter2)
        #             geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
        #             geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
        #         if self.small_space != 0:
        #             # small distance between centre of gratings in nm
        #             # calc complementary large distance, which is added to top & bottom
        #             large_d_on_2 = (self.period - self.diameter1/2 - self.diameter2/2 - self.small_space)/2
        #             posx1 = large_d_on_2 + self.diameter1/2
        #             posx2 = large_d_on_2 + self.diameter2/2
        #             posx3 = large_d_on_2 + self.diameter1 + ((self.small_space - self.diameter1/2 - self.diameter2/2)/2)
        #             geo = geo.replace('posx1 = hy/4;', "posx1 = %f/d_in_nm;" % posx1)
        #             geo = geo.replace('posx2 = hy/4;', "posx2 = %f/d_in_nm;" % posx2)
        #             geo = geo.replace('posx3 = hy/2;', "posx3 = %f/d_in_nm;" % posx3)
        #         # if supercell > 1:
        #         #     geo = geo.replace('a2 = 0;', "a2 = %i;" % self.diameter2)
        #         #     geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
        #         # if supercell > 2:
        #         #     geo = geo.replace('a3 = 0;', "a3 = %i;" % self.diameter3)
        #         #     geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)


        #         open(msh_location + msh_name + '.geo', "w").write(geo)
        #         EMUstack.conv_gmsh(msh_location+msh_name)
        #         # gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.msh'
        #         # gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.geo'
        #         # os.system(gmsh_cmd)

        else:
            raise ValueError("Must be simulating either a '1D_array' or a '2D_array'.")

    def calc_modes(self, light, **args):
        """ Run a simulation to find the NanoStruct's modes.

            Args:
                light  (Light instance): Represents incident light.

                args  (dict): Options to pass to :Simmo.calc_modes:.

            Returns:
                :Simmo: object
        """
        simmo = Simmo(self, light)

        simmo.calc_modes(**args)
        return simmo






class ThinFilm(object):
    """ Represents an unstructured homogeneous film.

        Args:
            period  (float): Artificial period imposed on homogeneous film \
                to give consistently defined plane waves in terms of \
                diffraction orders of structured layers.

        Keyword Args:
            period_y  (float): The period of the unit cell in the y-direction.\
                If None, period_y = period.

            height_nm  (float): The thickness of the layer in nm or 'semi_inf'\
                for a semi-infinte layer.

            num_pw_per_pol  (int): The number of plane waves per polarisation.

            world_1d  (bool): Does the rest of the stack have exclusively 1D \
                periodic structures and homogeneous layers? \
                If True we use the set of 1D diffraction order PWs.

            material  : A :Material: instance specifying the n of \
                the layer and related methods.

            loss  (bool): If False sets Im(n) = 0, if True leaves n as is.
    """
    def __init__(self, period, period_y=None, height_nm=1.0, num_pw_per_pol=0,
                 world_1d=False, material=materials.Material(3.0 + 0.001),
                 loss=True):
        self.period = float(period)
        if period_y is None:
            self.period_y = float(period)
        else:
            self.period_y = float(period_y)
        self.world_1d = world_1d
        self.height_nm = height_nm
        self.num_pw_per_pol = num_pw_per_pol
        self.material = material
        self.loss = loss

    def calc_modes(self, light):
        """ Run a simulation to find the ThinFilm's modes.

            Args:
                light  (Light instance): Represents incident light.

                args  (dict): Options to pass to :Anallo.calc_modes:.

            Returns:
                :Anallo: object
        """
        an = Anallo(self, light)
        an.calc_modes()
        return an






class Light(object):
    """ Represents the light incident on structure.

        Incident angles may either be specified by `k_parallel` or by
        incident angles `theta` and `phi`, together with the refractive
        index `n_inc` of the incident medium.

        `wl_nm` and `k_pll` are both in unnormalised units.

        At normal incidence and TE polarisation the E-field is aligned
        with the y-axis.

        At normal incidence some plane waves and Bloch modes become degenerate.
        This causes problems for the FEM solver and the ordering of the plane
        waves. To avoid this a small (1e-5) theta and phi are introduced.

        Args:

            wl_nm  (float): Wavelength, in nanometers.

        Keyword Args:
            max_order_PWs  (int): Maximum plane wave order to include.

            k_parallel  (tuple): The wave vector components (k_x, k_y) \
                parallel to the interface planes. Units of nm^-1.

            theta  (float): Polar angle of incidence in degrees.

            phi  (float): Azimuthal angle of incidence in degrees \
                measured from x-axis.
    """
    def __init__(self, wl_nm, max_order_PWs=2, k_parallel=None,
                 theta=None, phi=None, n_inc=1.):
        if np.imag(wl_nm) != 0:
            self.wl_nm = complex(wl_nm)
            print("Warning: using a complex wavelength. EMUstack can \n\
                only handle these for uniform films using 0 pw_orders.")
        else:
            self.wl_nm = float(np.real(wl_nm))
        self._air_anallos = {}
        self.max_order_PWs = max_order_PWs

        if None == theta and None == k_parallel:
            raise ValueError("Specify incident angle either by \n\
            k_parallel OR by theta, phi and n_inc.")

        if None == theta:
            self.k_pll = np.array(k_parallel, dtype='float64')
            # Check that not aligned with either x or y axis.
            if np.abs(self.k_pll[0]) == 0 or np.abs(self.k_pll[1]) == 0:
                print("Warning: a component of k_parallel is exactly zero, \n\
                this can lead to degeneracies and errors.")
        else:
            # Check for inconsistent input
            if None != k_parallel or phi == None:
                raise ValueError("Specify incident angle either by \n\
            k_parallel OR by theta, phi and n_inc.")
            # Avoid the degeneracies that occur at normal incidence
            # (FEM does not deal well with them)
            if abs(theta) < 1e-5: theta += 1e-5
            if abs(phi) < 1e-5: phi += 1e-5
            # Calculate k_parallel from incident angles
            k = 2 * np.pi * np.real(n_inc) / self.wl_nm
            theta *= np.pi / 180
            phi *= np.pi / 180
            self.k_pll = k*np.sin(theta) * np.array(
                        [np.cos(phi), np.sin(phi)], dtype='float64')



    def _air_ref(self, period, period_y, world_1d):
        """ Return an :Anallo: corresponding to this :Light: in free space.

            The :Anallo: will have len(anallo.k_z) == 2 * num_pw.

            Args:
                period  (float): period imposed on homogeneous film.

                period_y  (float): period imposed on homogeneous film \
                    along y-axis.

                world_1d  (bool): Specify whether to use 1D or 2D \
                    diffraction orders.
        """

        if (period) in self._air_anallos:
            return self._air_anallos[(period)]
        else:
            air = ThinFilm(period=period, period_y=period_y,
                           material=materials.Air, world_1d=world_1d)
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
    fmt_string = string.replace('.', '_')
    return fmt_string


def calculate_ff(inc_shape, d, dy, a1, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0,
                 a8=0, a9=0, a10=0, a11=0, a12=0, a13=0, a14=0, a15=0, a16=0,
                 el1=0):
    """ Calculate the fill fraction of the inclusions.

        Args:
            inc_shape  (str): shape of the inclusions.

            d  (float): period of structure, in same units as a1-16.

            dy  (float): period of structure along y-axis, in same units as a1-16.

            a1  (float): diameter of inclusion 1, in same units as d.

        Keyword Args:
            a2-16  (float): diameters of further inclusions.

            el1  (float): ellipticity of inclusion 1.
    """

    if inc_shape == 'circle' or inc_shape == 'ellipse':
        ff = np.pi*((a1/2)**2*np.sqrt(1-el1) + (a2/2)**2 + (a3/2)**2 +
                    (a4/2)**2 + (a5/2)**2 + (a6/2)**2 + (a7/2)**2 + (a8/2)**2 +
                    (a9/2)**2 + (a10/2)**2 + (a11/2)**2 + (a12/2)**2 + (a13/2)**2 +
                    (a14/2)**2 + (a15/2)**2 + (a16/2)**2)/(d*dy)
    elif inc_shape == 'square':
        ff = ((a1)**2 + (a2)**2 + (a3)**2 + (a4)**2 + (a5)**2 + (a6)**2 +
              (a7)**2 + (a8)**2 + (a9)**2 + (a10)**2 + (a11)**2 + (a12)**2 +
              (a13)**2 + (a14)**2 + (a15)**2 + (a16)**2)/(d*dy)
    elif inc_shape == 'dimer':
        ff = np.pi*((a1/2.0)**2+(a2/2.0)**2)/(d*dy)
    else:
        ff = 0.0
    return ff
