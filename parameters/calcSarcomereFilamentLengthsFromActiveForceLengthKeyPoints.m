function [activeForceLengthKeyPoints, ...
          halfMyosinBareLength, ...
          halfMyosinLength,...
          zLineLength,...
          actinLength] = ...
            calcSarcomereFilamentLengthsFromActiveForceLengthKeyPoints(animalId)
%%
%
% This function uses the method from Rassier et al. and Higuichi to 
% infer the length of the actin and myosin filament lengths from keypoints
% taken from the active force length relation of cat soleus muscle, human
% skeletal muscle, rabbit skeletal muscle and frog skeletal muscle
%
%  Rassier DE, MacIntosh BR, Herzog W. Length dependence of active force 
%  production in skeletal muscle. Journal of applied physiology. 
%  1999 May 1;86(5):1445-57.
%
%  Higuchi H, Yanagida T, Goldman YE. Compliance of thin filaments in skinned 
%  fibers of rabbit skeletal muscle. Biophysical journal. 1995 Sep 1;69(3):1000-10.
% 
% @param animalType
%   0. cat skeletal muscle
%   1. human skeletal muscle
%   2. frog skeletal muscle
%   3. rabbit skeletal muscle
%
% @return activeForceLengthKeyPoints: 
%     key points on the active force length curve in units of micro meters.
%
% @return halfMyosinBareLength:
%     half of the length of the bare patch on myosin in units of micro meters
%
% @return halfMyosinLength:
%     half of the total length of myosin in units of micro meters
%
% @return zLineLength:
%     half of the total thickness of the z-line in micro meters
%
% @return actinLength:
%     the length of one actin filament in micro meters
%%



halfMyosinBareLength = nan;
halfMyosinLength     = nan;                    
zLineLength          = nan;                    
actinLength          = nan;  

switch animalId

  case 1
    %From Rassier et al. Fig 3 & text
    geoCat   = [1.27, 1.7, 2.34, 2.51, 3.94];
    activeForceLengthKeyPoints = geoCat;  

    halfMyosinBareLength =     ( activeForceLengthKeyPoints(1,4) ...
                               - activeForceLengthKeyPoints(1,3) )*0.25;
    
    halfMyosinLength      = 0.5*(activeForceLengthKeyPoints(1,5) ...
                                -activeForceLengthKeyPoints(1,3));                    
    zLineLength           = 0.05;                    
    actinLength           = 0.5*( (activeForceLengthKeyPoints(1,5))...
                                    -2*halfMyosinLength ...
                                    -2*zLineLength);      

    here=1;
  case 2
    %From Rassier et al. Fig 3 & text    
    geoHuman = [1.27, 1.7, 2.64, 2.81, 4.24];    
    activeForceLengthKeyPoints = geoHuman;    

    halfMyosinBareLength =     ( activeForceLengthKeyPoints(1,4) ...
                               - activeForceLengthKeyPoints(1,3) )*0.25;
    
    halfMyosinLength      = 0.5*(activeForceLengthKeyPoints(1,5) ...
                                -activeForceLengthKeyPoints(1,3));                    
    zLineLength           = 0.05;                    
    actinLength           = 0.5*( (activeForceLengthKeyPoints(1,5))...
                            -2*halfMyosinLength ...
                            -2*zLineLength); 
    here=1;
    
  case 3
    %From Rassier et al. Fig 3 & text    
    % 0.95 um actin filament length
    % 1.6 um myosin filament length
    % 0.1 um zline width
    % 0.2 um bare length

    actinLength       = 0.95;
    myosinLength      = 1.6;
    zLineLength       = 0.05;
    myosinBareLength  = 0.2;

    halfMyosinBareLength = myosinBareLength*0.25;
    halfMyosinLength     = myosinLength*0.5;


    lasc  = max(actinLength,myosinLength)    + 2*zLineLength;    
    loptA = 2*actinLength                    + 2*zLineLength;
    loptB = 2*actinLength  +myosinBareLength + 2*zLineLength;
    lmax  = 2*actinLength  +myosinLength     + 2*zLineLength;

   
    geoFrog = [1.27,lasc, loptA, loptB,  lmax];  
    activeForceLengthKeyPoints = geoFrog;  

    flag_debug=0;
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

  case 4
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
  otherwise
    assert(0,'animalName must be cat, human, frog, or rabbit');
end



here=1;


          
