% VFIELD_HOLISTIC
%
% Usage:
%   [vest,Pvest] = vfield_holistic(vlos,dvlos,xgmag,ygmag,K,x2,y2)
%
% Inputs:
%    vlos  = line-of-sight velocity (from radar data)
%    dvlos = error standard deviation on vlos (from fitter)
%    xgmag = coordinates corresponding to vlos
%  / ygmag
%    K     = direction vectors [kx,ky,kz] (length(vlos) x 3)
%    Nx,Ny = solution grid size
%
% Outputs:
%    vest  : Estimated velocity field (Nx x Ny x 3)
%    Pvest : Covariance matrix for vest (3 x 3)
%
% Thomas Butler
% 6 Aug 2009


function [vest,Pvest,x2,y2] = vfield_holistic(vlos,dvlos,xgmag,ygmag,K,Nx,Ny)

sx = 500;
sy = 500;
sz = 15;

x2 = linspace(min(xgmag), max(xgmag), Nx+1);
y2 = linspace(min(ygmag), max(ygmag), Ny+1);

% Generate the geometry matrix
[A,Pv] = generate_matrices;

% Initialize covariance matrices
Pe = diag(dvlos.^2);

% Bayesian estimator for pixel (x,y)
vest = Pv*A'*inv(A*Pv*A' + Pe) * vlos;
Pvest = Pv - Pv*A'*inv(A*Pv*A' + Pe)*A*Pv;

% Reshape vest
vest = reshape(vest,[Ny,Nx,3]);

% Reshape Pvest
Pvx  = Pvest(         1:Ny*Nx ,          1:Ny*Nx );
Pvy  = Pvest(  Ny*Nx+(1:Ny*Nx),   Ny*Nx+(1:Ny*Nx));
Pvz  = Pvest(2*Ny*Nx+(1:Ny*Nx), 2*Ny*Nx+(1:Ny*Nx));
Pvxy = Pvest(         1:Ny*Nx ,   Ny*Nx+(1:Ny*Nx));
Pvxz = Pvest(         1:Ny*Nx , 2*Ny*Nx+(1:Ny*Nx));
Pvyz = Pvest(  Ny*Nx+(1:Ny*Nx), 2*Ny*Nx+(1:Ny*Nx));

Pvest = zeros(Ny,Nx,3,3);
for y = 1:Ny,
    for x = 1:Nx,
        Pvest(y,x,1,1) = Pvx(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,1,2) = Pvxy(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,1,3) = Pvxz(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,2,1) = Pvxy(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,2,2) = Pvy(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,2,3) = Pvyz(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,3,1) = Pvxz(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,3,2) = Pvyz(Ny*(y-1)+x,Ny*(y-1)+x);
        Pvest(y,x,3,3) = Pvz(Ny*(y-1)+x,Ny*(y-1)+x);
    end; % for y
end; % for x

x2 = ( x2(1:end-1) + x2(2:end) ) / 2;
y2 = ( y2(1:end-1) + y2(2:end) ) / 2;


% Nested function, generates the geometry matrix (A) and a priori
% covariance matrix (Pv) pixel-by-pixel.
%
% Size of A is [Nr,3*Ny*Nx]:
%
%  A = [ Kx | Ky | Kz ]
%
% where Ki represents the i-component of the direction vector for the beams
% in pixel (x,y).  The row index of A identifies the beam # (1 .. length of vlos).
% The column index corresponds to a column-major ordering of the x,y
% coordinates of vest (e.g. for 4x4, 1 = (y=1,x=1) .. 4 = (y=4,x=1) .. 16 = (y=4,x=4)).
%
% Size of Pv is [3*Ny*Nx,3*Ny*Nx]:
% 
%  Pv = [Pv_x |      |      ]
%       [     | Pv_y |      ]
%       [     |      | Pv_z ]
% 
% where the sub-matrices are square, positive semidefinite (covariance)
% matrices for each of the 3 spatial components.
%
% Thus, the result of vest = A*vlos is a [3*Ny*Nx,1] vector which can be
% reshaped simpy into three planes: reshape(vest,[Ny,Nx,3])
%
function [A,Pv] = generate_matrices
    
    A  = zeros(length(vlos),3*Ny*Nx);
    
    for x = 1:Nx,
        
        for y = 1:Ny,
            
            vlos_idx = find(xgmag >= x2(x) & xgmag <= x2(x+1) & ygmag >= y2(y) & ygmag <= y2(y+1));
%             vlos_current = vlos(vlos_idx);
%             dvlos_current = dvlos(vlos_idx);

            % Generate A matrix
            A(vlos_idx,Ny*(y-1)+x)           = K(vlos_idx,1); % x component
            A(vlos_idx,Ny*(y-1)+x + Ny*Nx)   = K(vlos_idx,2); % y component
            A(vlos_idx,Ny*(y-1)+x + 2*Ny*Nx) = K(vlos_idx,3); % z component
            
            % Generate Pv matrix
            Pvx = sx^2*eye(Ny*Nx);
            Pvy = sy^2*eye(Ny*Nx);
            Pvz = sz^2*eye(Ny*Nx);
            
            Pvxy = zeros(Ny*Nx);
            Pvxz = Pvxy;
            Pvyz = Pvxy;
            
            Pv = [ Pvx , Pvxy, Pvxz;
                   Pvxy, Pvy , Pvyz;
                   Pvxz, Pvyz, Pvz  ];
            
        end; % for x

    end; % for y

end % function generate_matrices


end % function vfield_holistic


