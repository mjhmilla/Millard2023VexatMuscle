function titinSegmentLengths = ...
	calcTitinSegmentLengths(sarcomereLength, ...
	                        cZToPEVKp, cZToPEVKd,...
	                        halfMyosinLength, ...
                            optimalSarcomereLength,...
                            zLineToT12Length)
%%
% This function takes the linear coefficents stored in cZToPEVKp and cZToPEVKd
% (returned by the function fitLinearModelToTitinsegmentElongation) and uses
% this function, along with the fixed filament lengths, to evaluate the 
% length of the prox-Ig, PEVK, and distal-Ig segments at a given sarcomere 
% length
%
% @param sarcomereLength the length of the sarcomere in units of micrometers
% @param cZToPEVKp 
% @param cZToPEVKd
% @param halfMyosinLength the length of the half-myosin in units of micrometers
% @param optimalSarcomereLength the length of the sarcomere at full myosin
%                               actin filament overlap.
% @return titnSegmentLengths: a struct with the fields:
%
%        (units of micrometers)
%        .lT12        :Z-line to T12 epitope 
%        .lIgp        :length of the prox-Ig segment
%        .lPevk       :length of the PEVK segment
%        .lIgdFree    :length of the distal-Ig segment that is not attached to myosin
%        .lIgdFixed   :length of the distal-Ig segment that is attahced to myosin
%        .lIgdTotal   :total length of the distal-Ig segment
%
%        (normalized by the optimal sarcomere length)
%        .lT12Norm
%        .lIgpNorm
%        .lPevkNorm
%        .lIgdFreeNorm
%        .lIgdFixedNorm
%        .lIgdTotalNorm
%%

titinSegmentLengths = ...
struct('lT12',nan, 'lIgp',nan, 'lPevk',nan, ...
	   'lIgdFree',nan,'lIgdFixed',nan,'lIgdTotal',nan,...
	   'lT12Norm',nan, 'lIgpNorm',nan, 'lPevkNorm', nan, ...
	   'lIgdFreeNorm',nan,'lIgdFixedNorm',nan,'lIgdTotalNorm',nan);



titinSegmentLengths.lT12       = zLineToT12Length;

titinSegmentLengths.lIgp       = [sarcomereLength, 1]*cZToPEVKp ...
                                 - titinSegmentLengths.lT12;

titinSegmentLengths.lPevk      = [sarcomereLength, 1]*cZToPEVKd ...
                                -[sarcomereLength, 1]*cZToPEVKp;

titinSegmentLengths.lIgdTotal  = 0.5*sarcomereLength ...
									-(  titinSegmentLengths.lT12...
									   +titinSegmentLengths.lIgp...
									   +titinSegmentLengths.lPevk);

titinSegmentLengths.lIgdFixed  = min(halfMyosinLength,...
	                                 titinSegmentLengths.lIgdTotal);

titinSegmentLengths.lIgdFree   = max(  titinSegmentLengths.lIgdTotal... 
  	                                  -titinSegmentLengths.lIgdFixed,0);

halfSarcomereLength = titinSegmentLengths.lT12 ...
                     +titinSegmentLengths.lIgp...
                     +titinSegmentLengths.lPevk...
                     +titinSegmentLengths.lIgdTotal;

assert(abs(halfSarcomereLength*2-sarcomereLength) < 1e-3);

titinSegmentLengths.lT12Norm       = ...
					titinSegmentLengths.lT12      / optimalSarcomereLength;

titinSegmentLengths.lIgpNorm       = ...
					titinSegmentLengths.lIgp      / optimalSarcomereLength;

titinSegmentLengths.lPevkNorm      = ...
					titinSegmentLengths.lPevk     / optimalSarcomereLength;

titinSegmentLengths.lIgdTotalNorm  = ...
					titinSegmentLengths.lIgdTotal / optimalSarcomereLength;

titinSegmentLengths.lIgdFixedNorm  = ...
					titinSegmentLengths.lIgdFixed / optimalSarcomereLength;

titinSegmentLengths.lIgdFreeNorm   = ...
					titinSegmentLengths.lIgdFree  / optimalSarcomereLength;

halfSarcomereLengthNorm = titinSegmentLengths.lT12Norm ...
                     +titinSegmentLengths.lIgpNorm...
                     +titinSegmentLengths.lPevkNorm...
                     +titinSegmentLengths.lIgdTotalNorm;

assert(abs(halfSarcomereLengthNorm*2 ...
      -(sarcomereLength/optimalSarcomereLength)) < 1e-3);
