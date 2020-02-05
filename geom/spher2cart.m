function [x,y,z]=spher2cart(r,theta,phi)
  z=r.*cos(theta);
  x=r.*sin(theta).*cos(phi);
  y=r.*sin(theta).*sin(phi);
end