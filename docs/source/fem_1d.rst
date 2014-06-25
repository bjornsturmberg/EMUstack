1D FEM Mode Solver
====================

1D Mesh
--------

1D FEM mesh are created by the python subroutine objects.make_mesh() and passed directly into the fortran routine 'py_calc_modes.f'.
The only parameter that influences this process is **'lc_bkg'**, where **1 / lc_bkg** is the number of FEM elements that the unit cell is divided into.

For a single inclusion the mesh is simply::

    |                   period                    |
        
    |--------------------| 1 |--------------------|

where the inclusion has 'diameter1' as is made of material 'inclusion_a'.

For a grating with 2 inclusions in the unit cell the spacing between the surfaces of the inclusions is set with the **'small_space'** parameter.::

    |                   period                    |
      | small_space  |      
    |2|----------------|   1   |--------------| 2 |

Inclusion1 will always be centered and of material 'inclusion_a', while all higher order inclusions are made of material 'inclusion_b'.

For unit cells that contain 3 or more inclusions there are 2 implemented spacing options.
By **default 'edge_spacing = False'** and the centers of all inclusions are equally spaced, with inclusion1 centered in the middle of the unit cell. ::

    |       |   equally    |  seperated   |       | 
    |--|    2    |------|  1  |---------| 3 |-----|

The alternative is to space the inclusions with equal distances between their surfaces. This is selected with the keyword argument **'edge_spacing = True'**::

    |  |       |equally |     |seperat-|   | ed   | 
    |--|   2   |--------|  1  |--------| 3 |------|

EMUstack can at present create mesh with up to 6 inclusions.
It is straightforward to extended this.
