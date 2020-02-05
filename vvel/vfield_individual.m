% VFIELD_INDIVIDUAL
%
% Usage:
%   [vest,Pvest] = vfield_individual(vlos,dvlos,xgmag,ygmag,K,Nx,Ny)
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


function [vest,Pvest,x2,y2] = vfield_individual(vlos,dvlos,xgmag,ygmag,K,Nx,Ny)

sx = 500;
sy = 500;
sz = 15;

x2 = linspace(min(xgmag), max(xgmag), Nx+2);
y2 = linspace(min(ygmag), max(ygmag), Ny+2);
    
vest = zeros(Ny,Nx,3);
Pvest = zeros(Ny,Nx,3,3);

for x = 1:Nx,
    for y = 1:Ny,

        % Grab the measurements corresponding to this "pixel"
        vlos_idx = find(xgmag >= x2(x) & xgmag <= x2(x+2) & ygmag >= y2(y) & ygmag <= y2(y+2));
        vlos_current = vlos(vlos_idx);
        dvlos_current = dvlos(vlos_idx);

        % Generate the geometry matrix
        A = K(vlos_idx,1:3);
        
        % Initialize covariance matrices
        Pe = diag(dvlos_current.^2);
        Pv = diag([sx, sy, sy].^2);

        % Bayesian estimator for pixel (x,y)
        vest(y,x,:) = Pv*A'*pinv(A*Pv*A' + Pe) * vlos_current;
        Pvest(y,x,:,:) = Pv - Pv*A'*pinv(A*Pv*A' + Pe)*A*Pv;

    end; % for y
end; % for x

x2 = x2(2:end-1);
y2 = y2(2:end-1);