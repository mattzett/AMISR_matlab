function duw=Warp_Allsky(fn_i,z0,lims)

% Map allsky image (fn_i) from az/el coordinates to x,y coordinates, where 
% x,y are horizontal cartesian coordinates in km at a specified altitude
% given by z0. A flat Earth is assumed. The image is flipped so that 
% 'contour' or 'quiver' can be correctly used to overlay fitted parameters.

% INPUTS:
% fn_i   Image to unwarp
% z0     Assumed altitude of the features (~120 km for aurora)
% lims   Limits for meshgrid [xlo,xhi,ylo,yhi]

% read image
d=fitsread(fn_i);

% Read el/az map for allsky images
path='/shared/classes/2008-2009/Fall2008/data/Allsky_Starfield_Fits/'; % Path on laptop
% path='~butler/classes/projects/AMISR-experiments/data/Allsky_Starfield_Fits/'; % Path on transcar server
fn_az='PKR_DASC_AZ_Nov2007.FITS';
az_i=fitsread([path fn_az]);
fn_el='PKR_DASC_EL_Nov2007.FITS';
el_i=fitsread([path fn_el]);

% NB: The following is done for FITS files of resolution 512x512 when 
%     the starfield fits az_i and el_i are 1024x1024.
if size(d) == [512,512],
    az_i = az_i(1:2:end,1:2:end);
    el_i = el_i(1:2:end,1:2:end);
end;

az_i = az_i * pi/180;
el_i = el_i * pi/180;

% flip image to accommodate apparent discrepancy between how matlab 
% implements 'imagesc' and 'contour' routines.
% df=fliplr(d);
% az_if=fliplr(az_i);
% el_if=fliplr(el_i);

% Create x and y pixel positions for altitude z0.
el_i(el_i<=10*pi/180)=10*pi/180;  %limit to elevations > 10 degrees
r=z0./sin(el_i);
x=r.*cos(el_i).*sin(az_i);
y=r.*cos(el_i).*cos(az_i);
fac=0.35;   %factor by which to reduce number of pixels to map.
xs=imresize(x,fac,'nearest');
ys=imresize(y,fac,'nearest');
ds=imresize(d,fac,'nearest');
% size(xs)
% size(ys)
% size(ds)

% Create data array over the same space
[xp,yp]=meshgrid(lims(1):1:lims(2),lims(3):1:lims(4));

% Interpolate onto regular xp,yp grid
duw=griddata(xs,ys,ds,xp,yp);

return
