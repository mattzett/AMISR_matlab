function [theta,phi]=azel2thetaphi(az,el)
  az2=mod(az(:)+360,360);
  phi=mod(90-az2+360,360);
  phi=phi*pi/180;
  
  theta=pi/2-el(:)*pi/180;
end