# SondePull
The code pulls radiosonde files from web archives and reads in variables for plotting and processing. It is a starter code for modification to meet a users needs.

## Authors
Dr. Elizabeth Smith, NOAA NSSL (elizabeth.smith@noaa.gov)

## Background
This algorthim was developed from some examples shared by Drs. Jeremy Gibbs and Kim Hoogewind of NSSL & OU-CIMMS.

It currently uses the Wyoming archive, but could be adapted to work with the IGRA2 arhcive as well! 


At present it reads in sounding data and uses the various methods (Mx) to detect PBL height defined in Seidel et al (2010): _Estimating climatological planetary boundary layer heights from radiosonde observations: Comparison of methods and uncertainty analysis_. JGR Atmos., *115(D16)*.

- M1 - parcel - height theta_v = theta_v(sfc)
- M1_a - parcel_conig - height theta_v = theta_v(sfc) + .6K from Coniglio et al (2013) Wea. Forecasting
- M2 - max pt gradient
- M3 - min spec. hum gradient
- M4 - min RH gradient
- M5 - min refractivity gradient -- do not have this varible so skipping
- M6 - base of elevated temp inversion
- M7 - top of sfc based temp inversion
- PBLh = median(M1-M7)

## Dependent pacakages

- datetime
- matplotlib
- metpy (https://github.com/Unidata/MetPy); note this code uses version 1.0 in which the specific_humidity_from_dewpoint function arguments are (pressure, dewpoint). Older versions have those arguments reversed. Check your version with <metpy.__version__> and either update or reverse the arguements as needed.
- netcdf4 -- could bypass if commenting out the writeout portion of the code.
- metpy (https://github.com/Unidata/MetPy)
- numpy
- pandas
- scipy
- siphon (https://github.com/Unidata/siphon)
