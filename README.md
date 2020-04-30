# neutrontools
Various utilities to work with neutron beam data

The pyhton module needed for reading SDF files can be found in  PaNOSC-ViNYL/SDF

In the folder SDF first cd to ./C/ folder and type “make”

Change directory to ../utilities and type “./build -s” - it will create a python module which is needed to read data from sdf files

Classes are defined to load and save particle data.

Dependencies: 
-------------
SDF module;
openPMD-api module;
SimEx package

*****************

See test2.py for usage
