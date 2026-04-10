# Global Gridded Crop Model (GGCM) Emulators implementation

## Description
Here is a Python module that provides class 'Emulator' that can be used
to instantiate crop-specific emulators of gridded crop growth models.
The emulator can estimate climatic average crop yields globally
on 0.5 arc-degree grid using polynomial coefficients as provided by and
described in the paper by
[Franke et al., 2020](https://doi.org/10.5194/gmd-13-3995-2020).

The the emulator is based on [AgMERRA Climate Forcing Dataset for
Agricultural Modeling](https://data.giss.nasa.gov/impacts/agmipcf/agmerra/)
available from NASA and described in
[Ruane et al., 2015](https://doi.org/10.1016/j.agrformet.2014.09.016).

An example of using the Emulator is included.

## Installation
1. Clone the repository.
2. Download a NetCDF file with emulator polynomial coefficients e.g. for
   [LPJmL model, maize crop](https://zenodo.org/records/3592453/files/LPJmL_maize_ggcmi_phase2_emulator_A0.nc4)
   and store it in the `data-emu` folder.
3. Start the emulator with `python ggcm_emulator.py`

## Usage
You would need to download more files containing emulator polynomial
coefficients to run `python usage_example.py`.
These can be obtained from [Zenodo](https://zenodo.org/records/3592453)
as described by
[Franke et al., 2020](https://doi.org/10.5194/gmd-13-3995-2020).

For details and usage example, please see the file `usage_example.py`.

## Acknowledgements

This work is supported by the Horizon Europe Research and Innovation
programme project DIAMOND (grant-no. 101081179).

## License
GNU General Public License
