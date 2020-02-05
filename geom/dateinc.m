function [UTnew,dmynew]=dateinc(UT,dmy,dt)

day=dmy(1); month=dmy(2); year=dmy(3);

UTnew=UT+dt/3600;
if UTnew>=24
    UTnew=mod(UTnew,24);
    daynew=day+1;
    
    switch month
        case {4,6,9,11}
            if daynew>30
                daynew=1;
                monthnew=month+1;
                monthinc=1;
            else
                monthnew=month;
                monthinc=0;
            end
        case 2
            if ~mod(year,4)
                daylim=28;
            else
                daylim=29;
            end
            
            if daynew>daylim
                daynew=1;
                monthnew=month+1;
                monthinc=1;
            else
                monthnew=month;
                monthinc=0;
            end
        otherwise
            if daynew>31
                daynew=1;
                monthnew=month+1;
                monthinc=1;
            else
                monthnew=month;
                monthinc=0;
            end
    end
    
    if monthinc & monthnew>12
        monthnew=1;
        yearnew=year+1;
    else
        yearnew=year;
    end
    
    dmynew=[daynew,monthnew,yearnew];
else
    dmynew=dmy;
end

end