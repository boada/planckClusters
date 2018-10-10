# This is the main directory with all of the code to reduce/analyze the MOSAIC1/2/3 and NEWFIRM imaging.

Here is a basic layout:

- bin -- executables for the dust correction
- bpz-1.99.3 -- the code for the BPZ redshifts -- this is massive and had to work with
- confs -- sextractor, stiff and swarp configuration files
- legacy -- old code I didn't want to throw away
- LIB -- code to do the dust corrections. In general this shouldn't need to be modified. You might need to recompile it for your specific machine if things don't work.

Here are the main working parts of the pipeline.

- pipe_utils -- high level scripts that the greater pipeline needs to work. Things like cosmology calculations, dust corrections, I/O, and other extras are in there. Most of theses will be imported as 'from pipe_utils import XXX'
- plugins -- these contain all of the working scripts that will go together to build the pipeline itself. Currently these are imported into the pipeline class. That makes all of these 'metaclasses.' These functions WILL NOT work if they are called or imported into another script. They are designed to work inside of the pipeline class ONLY. The file names generally describe which parts of the pipeline the functions correspond to. I chose to break everything up from a monolithic file to make bug tracking easier and to make it easier to navigate around in the pipeline.
- scripts -- these are utility scripts that are designed to be run on their own. They do big tasks like merge all the outputs together or move data around on the disks. Most of these will only need to be run once. There are some convinence scripts that I built to make things faster. The reduce_ALL_async script runs the reduction pipeline in parallel to speed everything up.


# The main pipeline file is combcat_PROJECTED.
This calls/imports all of the functions it needs from the folders above. It controls the setup and general handling of the pipeline, but it doesn't actually do any of the work. That is all handled in the other files.

- add_catalogs -- This file adds extra information to the calibrated sextractor catalogs before it goes into BPZ. Things like PS1/SDSS mags and spec-z's if they are around.
- utils -- Utility scripts that handle I/O for the pipeline or changing of data here and there. In general, there aren't any science scripts in there. The functions just do little things. 
