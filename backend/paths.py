"""
    paths.py contais all the relevant paths for the operation of the
    code. When path are assigned their default values, EMUstack code
    can only be executed from 'EMUstack/subfolder'. An appropriate
    modification of the path values with absolute paths allows for
    code execution from an arbitrary folder.



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

# where the python and fortran backend code is placed
backend_path = '../backend'

# folder containing relevant data as the optical constants
data_path = '../backend/data/'

# where template.geo files for mesh creation are store
template_path = '../backend/fortran/msh/'

# where existing .msh and .mail files are stored and
# new .msh and .mail files are created
msh_path = '../backend/fortran/msh/'
