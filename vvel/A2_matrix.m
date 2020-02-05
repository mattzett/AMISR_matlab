function A = A2_matrix(PART,az,el,dec,dip)
% A2_matrix
%  Generate the "A" matrix for transforming an array of vx/vy pairs into
%  line-of-sight projections.
% 
%  The geometry information is contained in the coordinate grids (X,Y) and
%  radar beam positions (xr,yr,zr) passed as inputs. Also accepts the
%  declination and dip angles (dec,dip).
% 
%   A = A2_matrix(PART,az,el,dec,dip)

[NxNy] = max(PART);
Nr      = length(az);

A = zeros(Nr, 2*NxNy);

% Rotation matrix
Rgmag = [cos(dec),         -sin(dec);
         sin(dip)*sin(dec), cos(dec)*sin(dip)];

% Count the number of dense grid points in each cell.
WEIGHT = zeros(NxNy,1);
for k = 1:(NxNy),
    WEIGHT(k) = numel(PART(PART==k));
end;

NORM = sqrt((cos(el).*sin(az)).^2 + (cos(el).*sin(az)).^2 + sin(el).^2);
for i = 1:Nr,
    % k vector for the current beam
    direction_vectors = [cos(el(i))*sin(az(i));
                         cos(el(i))*cos(az(i))];
    % Rotate k
%     direction_vectors = Rgmag*direction_vectors;
    % Scale by number of points in the cell
    direction_vectors = direction_vectors;

    if PART(i) > 0,
        start_col = 2*PART(i) - 1;
        end_col   = 2*PART(i);

        A(i,start_col:end_col) = direction_vectors(:);
    end;
end;
