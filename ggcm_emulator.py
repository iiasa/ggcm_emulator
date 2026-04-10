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


"""Global Gridded Crop Model (GGCM) Emulator implementation

What's in this module:
    Class 'Emulator' that can be used to instantiate crop-specific
    emulators of gridded crop growth models. The emulator can estimate
    crop yields globally on 0.5 arc-degree grid using polynomial
    coefficients as provided by and described in Franke et al., 2020,
    https://doi.org/10.5194/gmd-13-3995-2020

Typical use:
    import ggcm_emulator as emu
    emu.setup()
    e = emu.Emulator('maize')
    [Yield, W_oob, T_oob] = e.get_yields(C, T, W, N)

Usage example:
    See 'Usage example' at the end of this file.

:Author: Nikolay Khabarov, khabarov@iiasa.ac.at
:Date: 2026-04-09
"""

import numpy as np
import xarray as xr
import netCDF4 as netcdf

data_dir = './data-emu'
crop_model = 'LPJmL'
T_agmerra = 0             # placeholder for historical temperatures
W_agmerra = 0             # placeholder for historical precipitation
_setup_ok = False         # flag indicating successful setup

def setup(data_dir='./data-emu', crop_model='LPJmL'):
    """Function to carry out basic setup of this module

    Args:
        data_dir (string): data files folder;
        crop_model (string): crop model name.

    Returns:
        Nothing.

    Before crop-specific instances of emulators can be created, these
    actions have to be carried out by calling this function:

        1. Set location of data files folder 'data_dir', which contains:
            (a) Temperature and precipitation files from AgMERRA dataset
                (https://data.giss.nasa.gov/impacts/agmipcf/agmerra/)
                that were averaged over entire time period e.g.:
                  `ncra AgMERRA_{1980..2010}_prate.nc4 file1.nc4`
                and then aggregated from 0.25 arc-degree (native
                AgMERRA resolution) to 0.5 arc-degree (emulator
                resolution) e.g.:
                  `cdo gridboxavg,2,2 file1.nc4 file2.nc4`
                and then the coordinate ranges adjusted e.g.:
                  `cdo sellonlatbox,-180,180,90,-90 file2.nc4 final.nc4`
                The resulting file names for temperature and
                precipitation are respectively:
                  'agmerra-tavg-avg-1980-2010-05deg-adjlon.nc4'
                and
                  'agmerra-prate-avg-1980-2010-05deg-adjlon.nc4'.

            (b) Files containing polynomial coefficients specific to
                crop and global gridded crop growth model e.g.:
                  'LPJmL_maize_ggcmi_phase2_emulator_A0.nc4'.
                These files can be downloaded from a repository
                referenced in Franke et al., 2020,
                https://doi.org/10.5194/gmd-13-3995-2020

        2. Set global gridded crop growth model name 'crop_model' for
           later loading of polynomial emulator coefficients from files
           described in 1.(b) above.

        3. Load two data files described in 1.(a) above.

    The 'data_dir' and 'crop_model' settings are stored in the module
    namespace and can be accessed by a user. Changing them is allowed
    with caution and would lead to loading different data files (e.g.
    using different crop model) when crop-specific emulators are
    instantiated.

    The full list of global gridded crop models that can be used as the
    'crop_model' argument and corresponding datasets can be found in
    Franke et al., 2020, https://doi.org/10.5194/gmd-13-3995-2020
    """

    global T_agmerra, W_agmerra, _setup_ok

    _setup_ok = False # setup is yet incomplete

    # We want to access module global variables and assign them the
    # values of respective namesake local variables (function
    # arguments). These values can be accessed by a user via module
    # namespace.

    globals()['data_dir'] = data_dir
    globals()['crop_model'] = crop_model

    # Historical temperature and precipitation file names

    file_agmerra_ann_avg_Temp = (
        data_dir
        + '/agmerra-tavg-avg-1980-2010-05deg-adjlon.nc4'
    )

    file_agmerra_ann_avg_Prec = (
        data_dir
        + '/agmerra-prate-avg-1980-2010-05deg-adjlon.nc4'
    )

    # Loading data from netcdf files

    nc = netcdf.Dataset(file_agmerra_ann_avg_Temp, 'r')
    T_agmerra = nc.variables['tavg'][0,:,:]
    # units: degrees Celsius
    # (there is only one time value in the file, no need to keep three
    # dimensions)

    nc = netcdf.Dataset(file_agmerra_ann_avg_Prec, 'r')
    W_agmerra = nc.variables['prate'][0,:,:] * 365.25
    # converting units: mm/day to mm/year
    # (there is only one time value in the file, no need to keep three
    # dimensions)

    # Ensure that there are no zeros in 'W_agmerra' as we'll divide by
    # these values. Precipitation of < 1 mm/year is practically equal
    # to 1 mm/year) {$def:ref$1$}
    W_agmerra[W_agmerra < 1] = 1

    _setup_ok = True


class Emulator:
    """Crop-specific yield emulator

    The emulator is using 'data_dir' and 'crop_model' settings described
    in the module setup() function.
    """

    def __init__(self, crop_name):
        """Constructs crop-specific yield emulator

        This function is storing 'crop_name' in the namesake attribute
        of the class instance and loading polynomial coefficients from a
        netcdf file, using the settings configured by earlier call
        to module setup() function.
        """

        global _setup_ok, data_dir, crop_model

        if not _setup_ok:
            raise Exception("Can't instantiate emulator: setup() has \
                to be called first.")

        self.crop_name = crop_name

        # Keeping the data location for possible user debugging needs
        # as it contains current values of global variables 'data_dir'
        # and 'crop_model'.
        self.file_emulator_coefficients = (
            data_dir
            + '/' + crop_model
            + '_' + crop_name
            + '_ggcmi_phase2_emulator_A0.nc4'
        )

        # Load and store polynomial coefficients
        nc = netcdf.Dataset(self.file_emulator_coefficients, 'r')
        self.K = nc.variables['K_rf'][:,:,:]


    def get_yields(self, Ca, Ta, Wa, Na):
        """Estimates yield for a crop using polynomial emulator

        Args:
            Ca (int or float): CO2 concentration, global annual average,
                               in ppm, e.g.: 360.0

            Ta (numpy.array):  Temperature in a 0.5 degree cell, annual
                               average over e.g. 5 years, in degrees
                               Celsius, e.g.: 10.0

            Wa (numpy.array):  Precipitation in a 0.5 degree cell,
                               annual average over e.g. 5 years, in
                               mm/year, e.g.: 400.0

            Na (int or float): Fertilizer elemental N application,
                               kg/ha/year. Uniform for any location;
                               used e.g. for estimating yield potential,
                               e.g.: 200.0

        Each of the 'Ta' and 'Wa' arrays has:

            360 rows corresponding to pixel centroids latitudes
                from 89.75 (near North Pole; first row) to -89.75 (near
                South Pole; last row) at 0.5 arc-degree step;

            720 columns corresponding to pixel centroids longitudes
                from -179.75 (west of Greenwich; first column) to 179.75
                (east of Greenwich; last column) at 0.5 arc-degree step.

        Returns:
            List of three numpy arrays: [Yield, W_oob, T_oob] each in
            the same 360 rows and 720 columns format as this function
            arguments, where:

            Yield: estimated yield, in ton dry matter / year

            W_oob: water out of bounds value indicating (in mm/year)
                   by how much the respective 'Wa' element exceeds the
                   (historical-50%) .. (historical+30%) range valid for
                   application of the emulator;

            T_oob: temperature out of bounds value indicating (in
                   degrees Celsius) by how much the respective 'Ta'
                   element exceeded the (historical - 1) ..
                   (historical + 6) range valid for application of the
                   emulator;

            Negative values in W_oob and T_oob mean that the argument
            (respectively Wa or Ta) is below the lower bound of the
            range, positive values - above the upper bound.

        Remark:
            In case if Wa and/or Ta are out of bounds, the respective
            lower or upper bound is used to run the emulator to avoid
            potential issues with extrapolation.

            The Ca range of 360 .. 810 is respected in the same manner,
            but the out of bounds cases are not reported.
        """

        # Sanitize the arguments according to lower and upper bounds in
        # Franke et al., 2020 https://doi.org/10.5194/gmd-13-3995-2020

        C_san = min(max(360, Ca), 810)
        T_san = np.minimum(np.maximum(T_agmerra-1, Ta), T_agmerra+6)
        W_san = np.minimum(np.maximum(0.5*W_agmerra, Wa), 1.3*W_agmerra)
        N_san = min(max(10, Na), 200)

        ## "Out of bounds" arrays calculation
        T_oob = Ta - T_san
        W_oob = Wa - W_san

        # Shorter uniform naming of variables for using in polynomial;
        # for water input (precipitation) we need to carry out
        # conversion from absolute value (mm/year) to a multiplier W
        # required by emulator such that
        #   W_san = W * W_agmerra, where W in [0.5 .. 1.3],
        # that is
        #   W = W_san / W_agmerra.
        # Here we do not want to divide by zero, so we have already
        # prepared W_agmerra accordingly (see anchor in the code above
        # {$cite:ref$1$}).

        C = C_san
        T = T_san - T_agmerra   # we need difference from agmerra
        W = W_san / W_agmerra   # we need multiplier based on agmerra
        N = N_san
        K = self.K

        # Constructing polynomial. The index in K[] is shifted by -1
        # from that in the original Eq. (1) in Franke et al., 2020
        # (https://doi.org/10.5194/gmd-13-3995-2020), because
        # indexing in Python starts with 0, so that original K[1]
        # becomes K[0], original K[2] becomes K[1], etc.

        Yield = (
            K[0] + K[1]*C + K[2]*T + K[3]*W + K[4]*N + K[5]*C**2
          + K[6]*C*T + K[7]*C*W + K[8]*C*N + K[9]*T**2 + K[10]*T*W
          + K[11]*T*N + K[12]*W**2 + K[13]*W*N + K[14]*N**2
          + K[15]*C**3 + K[16]*C**2*T + K[17]*C**2*W + K[18]*C**2*N
          + K[19]*C*T**2 + K[20]*C*T*W + K[21]*C*T*N + K[22]*C*W**2
          + K[23]*C*W*N + K[24]*C*N**2 + K[25]*T**3 + K[26]*T**2*W
          + K[27]*T**2*N + K[28]*T*W**2 + K[29]*T*W*N + K[30]*T*N**2
          + K[31]*W**3 + K[32]*W**2*N + K[33]*W*N**2
        )

        # Replace all negative array elements with just zero, as yields
        # are always non-negative, yet polynomial form may occasionally
        # produce negative values.

        Yield[Yield < 0] = 0

        return [Yield, W_oob, T_oob]


def save05deg2nc(X, variable_name, units, file_name,
        description='Test output'):
    """Function for saving array to a NetCDF file (float32, compressed)

    Args:
        X (numpy.array): array of values that has:
            360 rows corresponding to pixel centroids latitudes
                from 89.75 (near North Pole; first row) to -89.75 (near
                South Pole; last row) at 0.5 arc-degree step;
            720 columns corresponding to pixel centroids longitudes
                from -179.75 (west of Greenwich; first column) to 179.75
                (east of Greenwich; last column) at 0.5 arc-degree step;

        variable_name (string): variable name to put into netcdf file;

        units (string): units of the variable to put into netcdf file;

        file_name (string): the name of netcdf file optionally including
            (full or relative) path to it;

        description (string): optional description to put into the file.

    Returns:
        Nothing.
    """

    # Define coordinates
    lat = np.arange(89.75, -90, -0.5)
    lon = np.arange(-179.75, 180, 0.5)

    # Create DataArray
    da = xr.DataArray(
        X,
        coords={'latitude': lat, 'longitude': lon},
        dims=['latitude', 'longitude'],
        name = variable_name
    )

    # Add metadata
    da.attrs['units'] = units
    da.attrs['description'] = description

    # Save to NetCDF, compression is {"zlib": True, "complevel": 9}
    parms = {variable_name:
                {"dtype": "float32", "zlib": True, "complevel": 9}}
    da.to_netcdf(file_name, encoding=parms)


########### Usage example. Runs if this module is executed. ############

if __name__ == "__main__":

    C = 360 # ppm
    T = np.ones((360, 720)) * 100 # degrees Celsius (high intentionally)
    W = np.ones((360, 720)) * 9000 # mm/year (high intentionally)
    N = 200 # kg/ha/year

    setup(data_dir = './data-emu', crop_model = 'LPJmL')
    e = Emulator('maize')

    [Yield, W_oob, T_oob] = e.get_yields(C, T, W, N)

    save05deg2nc(Yield, 'Yield', 't dm/ha/year', 'z-yield.nc4')
    save05deg2nc(W_oob, 'Water_out_of_bounds', 'mm/year', 'z-woob.nc4')
    save05deg2nc(T_oob, 'Temp_out_of_bounds', 'degrees C', 'z-toob.nc4')
