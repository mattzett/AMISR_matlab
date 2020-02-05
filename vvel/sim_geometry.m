function [X,Y,xr,yr,zr,PART] = sim_geometry(az,el,range,Nx,Ny,dec,dip)
% SIM_GEOMETRY
% Generate the geometry for a simulated atmospheric potential pattern.
%
% USAGE:
%  [X,Y,xr,yr,zr] = sim_geometry(az,el,range,Nx,Ny)
%  [X,Y,xr,yr,zr] = sim_geometry(az,el,range,Nx,Ny,dec,dip)
%
% INPUT:
%  az,el,range : Radar beam positions in spherical coordinates
%  Nx, Ny      : Number of elements in each direction for dense grid
%  dec, dip    : Declination and dip angles <optional>
%
% OUTPUT:
%  X, Y     : Mesh grid for potential pattern (dense) grid
%  xr,yr,zr : Vectors of radar beam directions (sparse)
%  PART     : Partition matrix, maps from (X,Y) to (xr,yr) <optional>

%----------------%
% Error checking %
%----------------%
if (nargin > 5),      % Input arguments include dec and dip angles?
    do_rotation = 1;
    if (nargin < 7),  % User forgot to include BOTH dec and dip
        error('SIM_GEOMETRY:inputErrorDecDip','Input argument must include both DEC and DIP, or neither.');
    end;
else
    do_rotation = 0;
end;

if (length(az) ~= length(el)),
    error('SIM_GEOMETRY:inputErrorAzEl','Az and El vectors must have the same length.');
end;

%--------------------------------%
% Convert spherical to cartesian %
%--------------------------------%
if (any(abs(el) > pi)) || any((abs(az) > 2*pi)),
    warning('SIM_GEOMETRY:anglesInRadians','Az > 2pi or El > pi.  Are angles in radians?')
end;
% Radar sampling points
xr = range.*cos(el).*sin(az);
yr = range.*cos(el).*cos(az);
zr = range.*sin(el);

if do_rotation,
    % Rotation matrix (geo -> gmag)
    Rgmag = [cos(dec),         -sin(dec),          0;
             sin(dip)*sin(dec), cos(dec)*sin(dip), cos(dip);
            -cos(dip)*sin(dec),-cos(dec)*cos(dip), sin(dip)];
    
    xyz = Rgmag*[xr(:)'; yr(:)'; zr(:)'];
    xr = xyz(1,:)';
    yr = xyz(2,:)';
    zr = xyz(3,:)';
end;

xmin = min(xr); ymin = min(yr);
xmax = max(xr); ymax = max(yr);

%----------------------------%
% Velocity field grid points %
%----------------------------%
x = linspace(xmin,xmax,Nx);
y = linspace(ymin,ymax,Ny);
[X,Y] = meshgrid(x,y);

if (nargout > 5),
%-------------------------------------------------------------------------%
% PARTITION the denser grid (X/Y) into cells corresponding to the         %
% center points of the sparser grid (xr/yr).                              %
%-------------------------------------------------------------------------%

    Nr = length(az);

    % Begin by computing the distance from the center of each sparse grid
    % element (xr/) to a point in the dense grid.
    STACK = zeros([Ny,Nx,Nr]);  % 3D matrix containing distances to each pt in the sparse grid (X/Y)
    for k = 1:Nr,
        STACK(:,:,k) = (X-xr(k)).^2 + (Y-yr(k)).^2;
    end;

    % Now for each element in PART, assign the index of the nearest element to
    % the sparse meshgrid.
    PART = zeros(Ny,Nx);
    for j = 1:Ny,
        for i = 1:Nx,
            PART(j,i) = find(STACK(j,i,:) == min(STACK(j,i,:)),1); % returns the desired index k
        end;
    end;

end;