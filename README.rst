Scripts used to produce the data and figures for the paper about the Arctic Ocean warming in CMIP5 models
=========================================================================================================

What's this?
------------

These are the scripts used to compute data and produce figures for the paper:

Burgard, C., and D. Notz (2017), Drivers of Arctic Ocean warming in CMIP5 models, Geophys. Res. Lett., 44, 4263–4271, https://doi.org/10.1002/2016GL072342.

Primary data
------------

The primary data was downloaded:
- from the ESGF website : https://esgf-data.dkrz.de/search/cmip5-dkrz/
- from the ceda ftp server : ftp://ftp.ceda.ac.uk/../badc/cmip5/data/cmip5/output1/
for the CMIP5 models listed in Table S1 in the Supplementary Information.

Computing data
--------------

The data was then processed using the Climate Data Operators cdo using following scripts:
- data_formatting.sh : "makes order" in the freshly downloaded files
- work_data_formatting.sh : prepares the data to use it further in python
- AOmask.sh: produces the mask for the Arctic Ocean and compute grid areas inside the Arctic Ocean
- different_regions.sh: produces the mask for the different individual regions
- meroceheatflux_MPI-ESM-LR.sh : computes the meridional oceanic heat flux directly from model output from MPI-ESM-LR (and not as a residual)
- understand_Hlat.sh: prepares data for further use in python that relates to explaining the changes in meridional oceanic heat flux
- AtlWat.sh: cuts out the upper 1000m to look at the Atlantic Water layer (especially in the Fram Strait)

Producing figures
-----------------

Final processing and visualization was done with python using following scripts:
- Figure1.py: creates Figure 1
- Figure2and4_rev.py : creates Figure 2 and 4
- Figure3_rev2.py: creates Figure 3

Signed: Clara Burgard, 26.05.2017
