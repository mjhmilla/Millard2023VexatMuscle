function [activeForceLengthKeyPoints, ...
          halfMyosinBareLength, ...
          halfMyosinLength,...
          zLineLength,...
          actinLength] = ...
            calcSarcomereFilamentLengthsFromActiveForceLengthKeyPoints(animalType)
%%
%
% This function uses the measurements and method from Rassier et al. to 
% infer the length of the actin and myosin filament lengths from keypoints
% taken from the active force length relation.
%
%  Rassier DE, MacIntosh BR, Herzog W. Length dependence of active force 
%  production in skeletal muscle. Journal of applied physiology. 
%  1999 May 1;86(5):1445-57.
%%

%The frog geometry has been blocked: the titin geometry below comes from
%a human soleus. Perhaps it is similar to cat skeletal muscel. However,
%there is no reason (that I) have to expect that it is similar to the titin
%molecule in frog skeletal muscle.
%
%geoFrog  = [1.27, 1.7,  2.0,  2.2, 3.6 ];

%Note that the length of the bare section is 1/4 of the length
%of the plateau. Why? The first half comes from the fact that 
%the section from geo(1,4) to geo(1,3) is a result of 2 actins moving
%across the flat section. We take half of that again because we're
%interested in half of the bare length.

geoCat   = [1.27, 1.7, 2.34, 2.51, 3.94];
geoHuman = [1.27, 1.7, 2.64, 2.81, 4.24];

switch animalType
  case 1
    activeForceLengthKeyPoints = geoCat;    
  case 2
    activeForceLengthKeyPoints = geoHuman;    
  otherwise
    assert(0, 'Error: flagFrog0Cat1Human1 incorrectly set');
end



halfMyosinBareLength =     ( activeForceLengthKeyPoints(1,4) ...
                           - activeForceLengthKeyPoints(1,3) )*0.25;

halfMyosinLength      = 0.5*(activeForceLengthKeyPoints(1,5) ...
                            -activeForceLengthKeyPoints(1,4)) ...
                            +halfMyosinBareLength;                    
zLineLength           = 0.05;                    
actinLength           = 0.5*( (activeForceLengthKeyPoints(1,5))...
                        -2*halfMyosinLength ...
                        -2*zLineLength);            
