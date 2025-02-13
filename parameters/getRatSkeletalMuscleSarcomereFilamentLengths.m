%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%
%%
function [activeForceLengthKeyPoints, ...
          halfMyosinBareLength, ...
          halfMyosinLength,...
          zLineLength,...
          actinLength] = getRatSkeletalMuscleSarcomereFilamentLengths()

% Extracted from the theoretical force-length relation in Figure 7 of
% Stephenson & Williams and Figure 1 of Tomalka:
%
%
%     Tomalka et al.    Stepenson & Williams    Zuurbier et a.
%       (um)            (um)
% l0:   1.3126                                  1.2939
% l1:   1.66            
% l2:   2.1288          2.516
% l3:   2.4013          2.6315
% l4:   3.9057          4.0
%
%           lo = 0.5(l2+l3) 
%       loT             loSW                    loZ
%       2.2651          2.5737                  2.38 (100-300 ms tetanic stim)
%
% l0N:  0.5795                                    
% l1N:  0.7329                          
% l2N:  0.9398          0.9776                 
% l3N:  1.0601          1.0225                  
% l4N:  1.7243          1.5542                
%
%
% * scaled from Tomalka et al to Stephenson & Williams using (loSW/loT)
%
% Actually Tomalka et al. reports filament lengths in the caption to Fig.1
%
%
% thick filament: 1.64 um  Zatsiorsky et al. 
% bare section  : 0.125    Gordon et al.
% thin filament : 1.13um   ter Keurs et al.
%
% Stephenson DG, Williams DA. Effects of sarcomere length on the force—pCa 
% relation in fast‐and slow‐twitch skinned muscle fibres from the rat. 
% The Journal of Physiology. 1982 Dec 1;333(1):637-53.
%
% Tomalka A, Heim M, Klotz A, Rode C, Siebert T. Ultrastructural and kinetic 
% evidence support that thick filaments slide through the Z-disc. Journal of 
% the Royal Society Interface. 2022 Dec 7;19(197):20220642.
%
% Gordon AM, Huxley AF, Julian FJ. 1966 The variation
% in isometric tension with sarcomere length in
% vertebrate muscle fibres. J. Physiol. 184, 170–192.
% (doi:10.1113/jphysiol.1966.sp007909).
%
% Zatsiorsky VM, Prilutsky BI. 2012 Biomechanics of
% skeletal muscles. Champaign, IL: Human Kinetics.
% (doi:10.5040/9781492595298)
%
% ter Keurs HE, Luff AR, Luff SE. 1984 Force —
% sarcomere-length relation and filament length in rat
% extensor digitorum muscle. Adv. Exp. Med. Biol. 170,
% 511–525. (doi:10.1007/978-1-4684-4703-3_44)
%
%Zuurbier CJ, Heslinga JW, Lee-de Groot MB, Van der Laarse WJ. Mean sarcomere 
%length-force relationship of rat muscle fibre bundles. Journal of Biomechanics. 
%1995 Jan 1;28(1):83-7.

geoRatSW = [1.4914,1.8862,2.516,2.63,4.0]; 
loSW   = 0.5*(geoRatSW(1,3)+geoRatSW(1,4));
geoRatSW = geoRatSW ./ loSW;


%From Higuchi et al. (cat)
% actinLength       = 1.12;
% myosinLength      = 1.63;
% zLineLength       = 0.07;
% myosinBareLength  = 0.16;

actinLength      = 1.13;    %ter Keurs et al.
myosinLength     = 1.63;    %Higuchi et al.
myosinBareLength = 0.125;   %Higuchi et al.
zLineLength      = 0.07;    %From Higuchi et al. (cat)
lzero            = 1.2939;  %From Zuurbier et al.

halfMyosinBareLength = myosinBareLength*0.5;
halfMyosinLength     = myosinLength*0.5;

lasc  = max(actinLength,myosinLength)    + 2*zLineLength;    
loptA = 2*actinLength  -myosinBareLength + 2*zLineLength;
loptB = 2*actinLength  +myosinBareLength + 2*zLineLength;
lmax  = 2*actinLength  +myosinLength     + 2*zLineLength;


geoRat = [lzero,lasc, loptA, loptB,  lmax]; 

activeForceLengthKeyPoints = geoRat;

flag_debug=1;
if(flag_debug==1)
    halfMyosinBareLengthTest = ( activeForceLengthKeyPoints(1,4) ...
                               - activeForceLengthKeyPoints(1,3) )*0.25;

    assert(abs(halfMyosinBareLengthTest-halfMyosinBareLength) < sqrt(eps));

    halfMyosinLengthTest = 0.5*( activeForceLengthKeyPoints(1,5) ...
                                -activeForceLengthKeyPoints(1,4)) ...
                                +halfMyosinBareLength;                    

    assert(abs(halfMyosinLengthTest-halfMyosinLength) < sqrt(eps));
                       
    actinLengthTest = 0.5*( (activeForceLengthKeyPoints(1,5))...
                            -2*halfMyosinLength ...
                            -2*zLineLength); 

    assert(abs(actinLengthTest-actinLength) < sqrt(eps));
    
    here=1;
end
here=1;