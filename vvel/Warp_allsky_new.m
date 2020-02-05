function [d_unwarp,lims]=Warp_Allsky(fn_img,z0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map Poker Flat allsky image, stored in FITS file fn_img, from az/el 
% to x,y coordinates, where x,y are horizontal cartesian coordinates mapped 
% to surface at altitude z0.  A flat Earth is assumed.  
%
% INPUTS:
% fn_i         string containing full path and filename of image
% z0           mapping altitude (km)
%
% OUTPUTS:
% duw          unwarped image
% lims         4-element vector containing bondaries of unwarped image
%              (xmin,xmax,ymin,ymax).  If the image is displayed using 
%              imagesc([lims(1),lims(2)],[lims(3),lims(4)],duw), then 
%              +y (downward direction) is North, and +x is East 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read el/az maps
fn_az='/shared/classes/projects/AMISR-experiments/data/Allsky_Starfield_Fits/PKR_DASC_AZ_Nov2007.FITS';
az_img=fitsread(fn_az);
fn_el='/shared/classes/projects/AMISR-experiments/data/Allsky_Starfield_Fits/PKR_DASC_EL_Nov2007.FITS';
el_img=fitsread(fn_el);

az_img = az_img * pi/180;
el_img = el_img * pi/180;

% read image
d=fitsread(fn_img);

% Compute cartesian coordinates from el,az,z0.  First limit el/az images 
% to consider only elevations greater some min_el (useful because griddata 
% is so slow).
min_el=30.*pi/180;          %minimum elevation to consider, in degrees
[goodx,goody]=find(el_img>=min_el);
el_img=el_img(min(goodx):max(goodx),min(goody):max(goody));
az_img=az_img(min(goodx):max(goodx),min(goody):max(goody));
d=d(min(goodx):max(goodx),min(goody):max(goody));
el_img(find(el_img<=min_el))=min_el;    %set elevations<min_el to a constant
r=z0./sin(el_img);
x=r.*cos(el_img).*sin(az_img);
y=r.*cos(el_img).*cos(az_img);
fac=.2;   %fraction by which to reduce number of pixels to map (griddata is really slow)
xs=imresize(x,fac,'nearest');
ys=imresize(y,fac,'nearest');
ds=imresize(d,fac,'nearest');

% Compute limits in output image
lims=[min(min(x)),max(max(x)),min(min(y)),max(max(y))];

% Create data array over the same space (1-km resolution)
[xp,yp]=meshgrid(lims(1):1:lims(2),lims(3):1:lims(4));

% Interpolate xs,ys onto regular xp,yp grid
d_unwarp=griddata(xs,ys,ds,xp,yp,'nearest');
el_d=atan(z0./sqrt(xp.^2+yp.^2));
d_unwarp(find(el_d<=min_el))=0;      %set pixel values below min_el to 0.

% imagesc([lims(1),lims(2)],[lims(3),lims(4)],d_unwarp);

return
