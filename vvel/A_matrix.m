function A = A_matrix(X,Y,xr,yr,zr,dec,dip)
% A_MATRIX
%  Generate the "A" matrix for transforming an array of vx/vy pairs into
%  line-of-sight projections.
%
% USAGE:
%  A = A_matrix(PART,az,el)
%  A = A_matrix(PART,az,el,dec,dip)
% 
% INPUT:
%  PART     : Partition matrix, mapping elements of X/Y to corresponding
%             positions in xr/yr
%  WEIGHT   : Number of X/Y elements in each direction; for weighting when
%             computing the mean.
%  az,el    : Radar beam directions
%  dec,dip  : Declination and dip angles <optional>
%
% OUTPUT:
%  A        : A matrix, mapping from coordinate meshgrids (X/Y, represented
%             by PART), to the radar beam positions (xr/yr, represented by
%             az & el).
%

%----------------%
% Error checking %
%----------------%
if (nargin > 5),      % Input arguments include dec and dip angles?
    do_rotation = 1;
    if (nargin < 7),  % User forgot to include BOTH dec and dip
        error('A_MATRIX:inputErrorDecDip','Input argument must include both DEC and DIP, or neither.');
    end;
else
    do_rotation = 0;
end;

if (length(xr) ~= length(yr)) || ...
   (length(yr) ~= length(zr)) || ...
   (length(zr) ~= length(xr)),
    error('A_MATRIX:inputErrorXrYrZr','xr/yr/zr vectors must have the same length.');
end;

%-----------------%
% Local variables %
%-----------------%
[Ny,Nx] = size(X);
Nr      = length(xr);

% If the following looks strange, it's because we're in East-from-North
% coordinates, and using elevation instead of zenith angle.
az = atan2(xr,yr);
el = atan2(zr,sqrt(xr.^2 + yr.^2));

% A = zeros(Nr, 2*Nx*Ny);
A = sparse(Nr, 2*Nx*Ny);

if do_rotation,
    % Rotation matrix (geo -> gmag)
    Rgmag = [cos(dec),         -sin(dec);
             sin(dip)*sin(dec), cos(dec)*sin(dip)];
else
    Rgmag = eye(2);
end;

for i = 1:Nr,
    % k vector for the current beam
    direction_vectors = [cos(el(i))*sin(az(i));
                         cos(el(i))*cos(az(i))];
    direction_vectors = Rgmag*direction_vectors;
    
    % Find the "pixel" index (X,Y) that is nearest (xr(i),yr(i))
    distances = (X-xr(i)).^2 + (Y-yr(i)).^2;
    j = find( distances == min(min(distances)) );
    
    % Assign direction vectors to appropriate columns of [A]_i.
    A(i,[2*j-1,2*j]) = direction_vectors(:);
end;

