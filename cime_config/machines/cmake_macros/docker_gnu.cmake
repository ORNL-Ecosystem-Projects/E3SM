# The local ELM Docker image installs the C and Fortran NetCDF libraries in a
# common prefix. These mixed-case variables are consumed by cprnc's CMake
# fallback when a NetCDF CMake package is not installed.
set(NetCDF_C_ROOT "$ENV{NETCDF_PATH}" CACHE PATH "NetCDF C prefix")
set(NetCDF_Fortran_ROOT "$ENV{NETCDF_PATH}" CACHE PATH "NetCDF Fortran prefix")
set(GENF90_PATH "${SRCROOT}/cime/CIME/non_py/externals/genf90" CACHE PATH "genf90 prefix")
