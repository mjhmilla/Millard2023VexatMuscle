function [activeForceLengthKeyPoints, ...
          halfMyosinBareLength, ...
          halfMyosinLength,...
          zLineLength,...
          actinLength] = getRabbitSkeletalMuscleSarcomereFilamentLengths()

%From Higuchi et al.
actinLength       = 1.12;
myosinLength      = 1.63;
zLineLength       = 0.07;
myosinBareLength  = 0.16;

halfMyosinBareLength = myosinBareLength*0.25;
halfMyosinLength     = myosinLength*0.5;

lasc  = max(actinLength,myosinLength)    + 2*zLineLength;    
loptA = 2*actinLength                    + 2*zLineLength;
loptB = 2*actinLength  +myosinBareLength + 2*zLineLength;
lmax  = 2*actinLength  +myosinLength     + 2*zLineLength;

geoRabbit = [1.27,lasc, loptA, loptB,  lmax]; 

activeForceLengthKeyPoints = geoRabbit;

flag_debug=1;
if(flag_debug==1)
    halfMyosinBareLength =     ( activeForceLengthKeyPoints(1,4) ...
                               - activeForceLengthKeyPoints(1,3) )*0.25;
    
    halfMyosinLength      = 0.5*(activeForceLengthKeyPoints(1,5) ...
                                -activeForceLengthKeyPoints(1,3));                    
                       
    actinLength           = 0.5*( (activeForceLengthKeyPoints(1,5))...
                            -2*halfMyosinLength ...
                            -2*zLineLength); 
end
here=1;