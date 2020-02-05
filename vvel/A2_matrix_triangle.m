function A = A2_matrix_triangle(VERTICES,az,el,dec,dip)
% A2_matrix_triangle
%  Generate the "A" matrix for transforming an array of vx/vy pairs into
%  line-of-sight projections.
% 
%  The geometry information is contained in the triangle VERTICES and
%  radar beam positions (az,el) passed as inputs. Also accepts the
%  declination and dip angles (dec,dip).
% 
%   A = A2_matrix_triangle(VERTICES,az,el,dec,dip)

% Error checking
if (size(VERTICES,2) ~= 3),
    error('VERTICES must be and Nx3 matrix of triangle vertices');
end;

N  = size(VERTICES,1);
Nr = length(az);

A = zeros(Nr, 2*N);

% k vector for the current beam
direction_vectors = [cos(el).*sin(az), cos(el).*cos(az)];

% NORM = sqrt((cos(el).*sin(az)).^2 + (cos(el).*sin(az)).^2 + sin(el).^2);
for i = 1:N,
    for j = VERTICES(i,:),
        A(j,(2*i-1):(2*i)) = direction_vectors(j,:);
    end;
end;
