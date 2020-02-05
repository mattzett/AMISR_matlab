function [Xenu,Yenu,Zenu] = Geodetic2ENU(GeodLat,GeodLong,GeodAlt,viewLat,viewLong)
    
% convert this to ECR
    a = 6378.1;
    b = 6356.8;
    e = sqrt(a^2-b^2)/a;
    N = a./( sqrt(1-e^2.*sind(GeodLat)));
    
    X = (N+GeodAlt).*cosd(GeodLat).*cosd(GeodLong);
    Y = (N+GeodAlt).*cosd(GeodLat).*sind(GeodLong);
    Z = ( N*(1-e^2)+GeodAlt ).*sind(GeodLat);
        
% find viewing geocentric coords from geodetic
    % viewing XYZ
    No = a./( sqrt(1-e^2.*sind(viewLat)));
    Xo = (No+.213).*cosd(viewLat).*cosd(viewLong);
    Yo = (No+.213).*cosd(viewLat).*sind(viewLong);
    Zo = (No*(1-e^2)+.213 ).*sind(viewLat);
    
% Center the x,y,z location based on viewing location
    Xctr = X - Xo;
    Yctr = Y - Yo;
    Zctr = Z - Zo;
    
    viewLat = 180/pi*atan2(Zo,sqrt(Xo.^2+Yo.^2));
    viewLong = 180/pi*atan2(Yo,Xo);
    
% convert from GEO to ENU 
    ECEF_to_ENU_RM = [-sind(viewLong)               , cosd(viewLong)                 , 0;
                     -sind(viewLat)*cosd(viewLong), -sind(viewLat)*sind(viewLong), cosd(viewLat);
                     cosd(viewLat)*cosd(viewLong) , cosd(viewLat)*sind(viewLong) , sind(viewLat)];
                 
    Xenu = ECEF_to_ENU_RM(1,1)*Xctr + ECEF_to_ENU_RM(1,2)*Yctr + ECEF_to_ENU_RM(1,3)*Zctr;
    Yenu = ECEF_to_ENU_RM(2,1)*Xctr + ECEF_to_ENU_RM(2,2)*Yctr + ECEF_to_ENU_RM(2,3)*Zctr;
    Zenu = ECEF_to_ENU_RM(3,1)*Xctr + ECEF_to_ENU_RM(3,2)*Yctr + ECEF_to_ENU_RM(3,3)*Zctr;
    
end