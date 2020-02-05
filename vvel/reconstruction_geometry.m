function [VERTICES,X2,Y2,Z2] = reconstruction_geometry(xr,yr,zr,type)
% RECONSTRUCTION_GEOMETRY
%  Generate the x,y,z coordinates of cell centers for the reconstruction.
%
% USAGE:
%  [VERTICES,X2,Y2,Z2] = reconstruction_geometry(type)
%
% INPUTS:
%  type     : either "triangle" or "square"
%
% OUTPUTS:
%  VERTICES : vertices of the partitioning
%  X2,Y2,Z2 : coordinates of the centroid of each cell
%
% Thomas Butler
% 14 April, 2009

if strcmp(type,'triangle'),
    
    VERTICES = delaunay(xr,yr);
    X2 = sum(xr(VERTICES),2)/3;
    Y2 = sum(yr(VERTICES),2)/3;
    Z2 = sum(zr(VERTICES),2)/3;

elseif strcmp(type,'square'),
    
    VERTICES = [ 1  2 22 25;
                 2  3 20 22;
                 3 16 17 20;
                16  4 14 17;
                25 22  9  8;
                22 20 10  9;
                20 17 11 10;
                17 14 12 11;
                 8  9 23 26;
                 9 10 21 23;
                10 11 18 21;
                11 12 15 18;
                26 23 24  5;
                23 21  6 24;
                21 18 19  6;
                18 15  7 19];
    X2 = sum(xr(VERTICES(:,[1,3])),2)/2;
    Y2 = sum(yr(VERTICES(:,[1,3])),2)/2;
    Z2 = sum(zr(VERTICES(:,[1,3])),2)/2;
    
else
    
    error('RECONSTRUCTION_GEOMETRY:InvalidType','type must be either ''triangle'' or ''square''');
    
end;