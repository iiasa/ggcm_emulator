#!/usr/bin/python

# Copyright (c) 2026 International Institute for Applied Systems
#                    Analysis (IIASA)
#
# Copyright (c) 2026 Nikolay Khabarov
#
# This is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This software is made available in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You can obtain a copy of the GNU General Public License at
# https://www.gnu.org/licenses/gpl-3.0.en.html or by writing to
# the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


"""This is an example of using Global Gridded Crop Model (GGCM) Emulator

:Author: Nikolay Khabarov, khabarov@iiasa.ac.at
:Date: 2026-04-09
"""

import ggcm_emulator as emu
import numpy as np


############### Define inputs for estimating crop yields ###############

# Global annual average CO2 concentration, ppm
C = 360

# Gridded annual average temperature, degrees Celsius.
# These high values imply using emulator's upper bound for temperature.
T = np.ones((360,720)) * 100

# Gridded annual precipitation/water, mm/year.
# These high values imply using emulator's upper bound for
# precipitation.
W = np.ones((360,720)) * 9000

# Fertilizer elemental N application, kg/ha/year (uniform for any
# location; used e.g. for estimating yield potential).
N = 200

# Each of the T and W arrays has:
#   360 rows corresponding to pixel centroids latitudes from 89.75
#       (near North Pole; first row) to -89.75 (near South Pole;
#       last row) at 0.5 arc-degree step;
#   720 columns corresponding to pixel centroids longitudes from -179.75
#       (west of Greenwich meridian; first column) to 179.75
#       (east of Greenwich meridian; last column) at 0.5 arc-degree
#       step.


######################## Emulator module setup #########################

emu.setup(data_dir='./data-emu', crop_model='LPJmL')


############### Running emulators for a list of crops ##################

crop_list = ['maize', 'rice', 'soy', 'spring_wheat', 'winter_wheat']
emu_list = []

# Create crop-specific emulators and add them to the list
for crop_name in crop_list:
    emu_list.append(emu.Emulator(crop_name))

for e in emu_list:
    # Apply emulator to estimate yields globally at 0.5 arc-degree
    [Yield, W_oob, T_oob] = e.get_yields(C, T, W, N)

    # Save the yields in a NetCDF file
    emu.save05deg2nc(
        Yield,

        variable_name='Yield_' + e.crop_name + '_' + emu.crop_model,

        units='t dm/ha/year',

        file_name='./z-test-yield_' + e.crop_name + '_' + emu.crop_model
            + '_C' +str(C)+ '_N' +str(N)+ '.nc4'
    )

    # The information on "out of bounds" values of T and W contained in
    # T_oob and W_oob is ignored here, but can be analyzed if needed.
