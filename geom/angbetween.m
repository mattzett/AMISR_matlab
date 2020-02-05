function ang=angbetween(az1,el1,az2,el2)

[theta1,phi1]=azel2thetaphi(az1,el1);
r1=1;
[theta2,phi2]=azel2thetaphi(az2,el2);
r2=1;

pos1=[r1*sin(theta1)*cos(phi1); r1*sin(theta1)*sin(phi1); r1*cos(theta1)];
pos2=[r2*sin(theta2)*cos(phi2); r2*sin(theta2)*sin(phi2); r2*cos(theta2)];

ang=acos(pos1'*pos2);

end