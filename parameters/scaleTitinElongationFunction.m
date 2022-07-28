function [lineZToPevkP, lineZToPevkD] = ...
    scaleTitinElongationFunction(...
            optimalSarcomereLength, lT12, halfMyosinLength,...
                numDomainsIgP, numResiduesPevk, numDomainsIgD, ...
            optimalSarcomereLengthRef, lT12Ref, halfMyosinLengthRef, ...
                numDomainsIgPRef, numResiduesPevkRef, numDomainsIgDRef,...
            maxIgDomainStrain_um,maxPevkResidueStrain_um,...
            lineZToPevkPRef, lineZToPevkDRef)
%%
% This function will scale the linear functions lineZToPevkPRef, lineZToPevkDRef
% that describe how the distance between the z-line and the proximal end of
% the PEVK segment (and distal end) varies with the length of the sarcomere.
% In brief, this function will adjust the proportion of the prox. Ig, PEVK,
% and distal Ig segments in proportion to the change in segment lengths relative
% to the refrence. Thus if the reference titin element has 
%
%   68 prox. Ig domains
%   2174 PEVK residues
%   22 distal Ig domains 
% 
% and the target has 
%
%   50 prox. Ig domains
%   800 PEVK residues
%   22 distal Ig domains 
%
% Then the PEVK segment will elongate in proportion to the ratio of
% [PEVK/(IgP+IgD)] / [PEVKRef/(IgPRef+IgDRef)]
%
% Similarly, the IgP segment will elongate in proportion to 
% [IgP/(IgP+IgD)] / [IgPRef/(IgPRef+IgDRef)]
%
%%

ratioPevkIg     = numResiduesPevk    / (numDomainsIgP+numDomainsIgD);
ratioPevkIgRef  = numResiduesPevkRef / (numDomainsIgPRef+numDomainsIgDRef);

ratioIgP        = numDomainsIgP / (numDomainsIgP+numDomainsIgD);
ratioIgPRef     = numDomainsIgPRef / (numDomainsIgPRef+numDomainsIgDRef);

ratioIgD        = numDomainsIgD / (numDomainsIgP+numDomainsIgD);
ratioIgDRef     = numDomainsIgDRef / (numDomainsIgPRef+numDomainsIgDRef);


%Construct linear functions to evaluate the length of the prox Ig, Pevk, and
%distal Ig segments


linePevkRef     = lineZToPevkDRef - lineZToPevkPRef;

%The lines of best fit to the data of Trombitas et al. end up suggesting
%that the amount of strain per Ig domain differs in the proximal and distal
%segments: 0.0028 / domain at lopt (prox.) vs 0.0038 at lopt (distal). 
%Which is strange, but somewhat expected due to the difficulty of the 
%experiments. Here we are first going to re evaluate these strain 
%functions such that the strain per Ig domain is uniform.
lineIgRefTotal = [0.5;0] - linePevkRef ...
                 -[0;halfMyosinLengthRef] - [0;lT12Ref];

%lineIgPRef      = lineZToPevkPRef - [0;lT12Ref];
lineIgPRef = lineIgRefTotal.*ratioIgPRef;
lineIgDRef = lineIgRefTotal.*ratioIgDRef;


%lineIgDRef     = [0.5;0] - lineZToPevkDRef - [0;halfMyosinLengthRef];
%lineIgDRef      = [0.5;0] - linePevkRef   - lineIgPRef ...
%                 -[0;halfMyosinLengthRef] - [0;lT12Ref];
%lineIgPRefTest = lineIgTotalRef.*ratioIgPRef;
%lineIgDRefTest = lineIgTotalRef.*ratioIgDRef;

loptIgPRef  = [optimalSarcomereLengthRef, 1]*lineIgPRef;
loptPevkRef = [optimalSarcomereLengthRef, 1]*linePevkRef;
loptIgDRef  = [optimalSarcomereLengthRef, 1]*lineIgDRef;
loptHalf = lT12Ref + loptIgPRef + loptPevkRef + loptIgDRef...
          + halfMyosinLengthRef;
assert( abs(loptHalf - 0.5*optimalSarcomereLengthRef) < 1e-3 );

%Check that both Ig segments have the same length per domain
loptIgPRefDomainStrain = loptIgPRef/numDomainsIgPRef;
loptIgDRefDomainStrain = loptIgDRef/numDomainsIgDRef;
loptPevkRefResidueStrain= loptPevkRef/numResiduesPevkRef;

assert(abs( (loptIgPRefDomainStrain)-(loptIgDRefDomainStrain)) < 1e-3 );

%Check that the domain and residue lengths are below the contour lengths
%maxIgDomainStrain_um = 25/1000;    %25 nm
%maxPevkResidueStrain_um = 0.38/1000; %0.38 nm
assert( loptIgPRefDomainStrain < maxIgDomainStrain_um);
assert( loptPevkRefResidueStrain < maxPevkResidueStrain_um);

%Scale the reference titin model to the target
linePevk    = linePevkRef.*(ratioPevkIg/ratioPevkIgRef);
lineIgTotal = [0.5;0] - linePevk - [0;halfMyosinLength] - [0;lT12];

lineIgP     = lineIgTotal.*(ratioIgP);
lineIgD     = lineIgTotal.*(ratioIgD);

lineZToPevkP = lineIgP      + [0;lT12];
lineZToPevkD = lineZToPevkP + linePevk;


loptIgP  = [optimalSarcomereLength, 1]*lineIgP;
loptPevk = [optimalSarcomereLength, 1]*linePevk;
loptIgD  = [optimalSarcomereLength, 1]*lineIgD;
loptHalf = lT12 + loptIgP + loptPevk + loptIgD + halfMyosinLength;
assert( abs(loptHalf - 0.5*optimalSarcomereLength) < 1e-3 );


%Check that both Ig segments have the same length per domain
loptIgPDomainStrain = loptIgP/numDomainsIgP;
loptIgDDomainStrain = loptIgD/numDomainsIgD;
loptPevkResidueStrain= loptPevk/numResiduesPevk;
assert(abs( (loptIgPDomainStrain)-(loptIgDDomainStrain)) < 1e-3 );

%Check that the domain and residue lengths are below the contour lengths
assert( loptIgPDomainStrain < maxIgDomainStrain_um);
assert( loptPevkResidueStrain < maxPevkResidueStrain_um);


here=1;


%  Amount of PEVK stretch will vary with ratio of PEVK to Ig
%     ( numResiduesPevk/(numDomainsIgP + numDomainsIgD) )
%   / ( numResiduesPevkRef/(numDomainsIgPRef+ numDomainsIgDRef) )
%
%Amount of IgP and IgD stretch will vary with
% (numDomainsIgP/numDomainsIgD) / (numDomainsIgPRef/numDomainsIgDRef)
