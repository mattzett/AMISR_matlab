# AMISR_matlab
Matlab scripts for reading and plotting AMISR data

## DESCRIPTIONS OF AMISR CODES

1 - amisr
The master code, calls and runs all following codes sequentially and creates main outputs. User must specify the file name, denoted "filelab", which will populate through the remainder of the code accordingly. The start and stop time steps ("itstart" and "itfin") can be used for unique time window specifications. 


2 - hdf5_extract
This function takes in PFISR HDF5 files, as found on amisr.com/somewhere..., and outputs the data in matlab-friendly form (_rawdata.mat). 

3 - plot_grid
Using the data loaded from amisr, plots a series of PNG images that display multiple dotted blue lines (each indicating a radar beam), black lines (indicating the magnetic field at the location of PFISR), and the source of the radar as a blue circle at the bottom of the plot. Rocket trajectories can be added to this plot using the trajectory flag and the resulting image is saved.

4 - grid_ISR
Converts the radar beam location from az and el coordinates to 3D cartesian wrt PFISR at the origin. Creates magnetic field lines too, spawning from ~150 km on each radar beam and is roated to the appropriate magnetic field line orientation. This grid information is saved (_fieldgrid.mat).

5 - makemovie_AMISR_vlos  
Using data from plot_grid and grid_ISR, creates plots of the usual plasma parameters in a configuration of 3-D slices, which are then used to create a movie of all individual images. The method used to make the movie depends on the operating system running matlab, a flag needs to be set.

6 - estimate_flowfield
Calculates the flow fields using the Butler 2010 method and saves the processed data (_vvelscans.mat). Creates a keogram of eastward and northward drift velocities, expressed as "distance north" from PFISR in km, against a horizontal axis of UT. 

7 - plot_ISRflows
Creates individual subplots of isti, isne, and fitted velocity magnitude at a specific altitude currently set to 300km, as well as the flow fields plotted using arrows in the left panel.


## DESCRIPTION OF OTHER PLOTTING SCRIPTS

1 - ButlerFits
Loads "_vvelscans.mat" and "_fieldgrid.mat" data, creates a keogram of eastward and northward drift velocities, expressed as "distance north" from PFISR in km, against a horizontal axis of UT. Utilizes Butler 2010 method.

2 - VvelFits
Loads vvel dataset from PFISR, creates a keogram of eastward and northward drift velocities, expressed as magnetic latitude "mlat." from PFISR, against a horizontal axis of UT. Utilizes standard PFISR method.

3 - plot_OneParamAllBeams
Intakes "_rawdata.mat, _fieldgrid.mat, _vvelscans.mat" and creates figure displaying all beams, isolating one parameter to observe. Parameters include "isne, isp, iste, isti, isvi". Uses AlphaData to make any extreme errors transparent in the final figure.

4 - plot_AllParamsOneBeam
Intakes "_rawdata.mat", creates two figures, one to display the parameters of "isne, isti, iste, isvi", and a second to display the corresponding errors of each parameter. Uses AlphaData to make any extreme errors transparent in the final figure.

