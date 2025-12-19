%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
% If you use this code in your work please cite the pre-print of this paper
% or the most recent peer-reviewed version of this paper:
%
%    Matthew Millard, David W. Franklin, Walter Herzog. 
%    A three filament mechanistic model of musculotendon force and impedance. 
%    bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 
%
%%

function [lineZToPevkP, lineZToPevkD] = ...
    scaleTitinElongationFunction2025(...
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
%Which is strange, but (perhaps?) somewhat expected due to the difficulty  
%of the experiments. Here we are first going to re evaluate these strain 
%functions such that the strain per Ig domain is uniform.
lineIgRefTotal = [0.5;0] - linePevkRef ...
                 -[0;halfMyosinLengthRef] - [0;lT12Ref];

%lineIgPRef      = lineZToPevkPRef - [0;lT12Ref];
lineIgPRef = lineIgRefTotal.*ratioIgPRef;
lineIgDRef = lineIgRefTotal.*ratioIgDRef;

lineIgPExact = (lineZToPevkPRef-[0;lT12Ref]);
lineIgDExact = ([0.5;0]-[0;halfMyosinLengthRef]-lineZToPevkDRef);

here=1;
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


%
nIg      = numDomainsIgP+numDomainsIgD;
nIgRef   = numDomainsIgPRef+numDomainsIgDRef;
nPevk    = numResiduesPevk;
nPevkRef = numResiduesPevkRef;

%
% As before, I'm lumping the Ig segment into a spring and the PEVK 
% segment into a spring. Our goal here is to evaluate the ratio
% of the length of the target PEVK segment to the reference PEVK 
% segment length.
%
% We start by noting that for two springs in series (kig1 and kp1) the
% length under force 1 is 
%
% l = f1/kig1 + f1/kp1
%
% For our target titin filament, we also will stretch it to a length
% l, but of course this will require an unknown force f2
%
% l = f2/kig2 + f2/kp2 
%
% Assuming that each Ig domain and PEVK residue in each filament has the
% same stiffness in both the reference and target domains, and that the
% stiffness is linear in the region of interest
%
% kig1 = (1/Ni1)*kig
%
% where kig1 is the stiffness of the entire proximal and distal Ig 
% segment (as if they were connected in series), Ni1 is the total number
% of Ig domains in the lumped segment, and kig is the stiffness of a single
% domain. Similarly for our target titin
%
% kig2 = (1/Ni2)*kig
%
% And also for the PEVK segments
%
% kp1 = (1/Np1)*kpevk
% kp2 = (1/Np2)*kpevk
%
% of the two titin filaments
%
% Finally we also need
%
% kig/kpevk.
%
% Under the assumption of constant stiffness
%
% kig = delta f / delta lig
%
% and
%
% kpevk = delta f/ delta lpevk
%
% Since delta f is unknown we cannot calculate kig nor kpevk.
% However, from Trombitas data, we can calculate the ratio:
%
% kig/kpevk = (1/delta lig)/(1/delta lpevk)
%
% since the delta f's cancel out. Since we know the stretch rate of the 
% Ig segments and the PEVK segment from Trombitas et al., and we also
% know the number of Ig domains and PEVK residues, we can evaluate
% kig/kpevk
%
igStretchRate   = lineIgRefTotal(1,1)*(1/nIgRef);
pevkStretchRate = linePevkRef(1,1)*(1/nPevkRef);
kig_kpevk       = pevkStretchRate/igStretchRate;

% Since the stiffness of the Ig segment is just the series stiffness of
% Nig domains in series we have
% 
% kig1 = kig/Ni1
% kig2 = kig/Ni2
% 
% where kig is the stiffness of an individual Ig domain.
% Similarly
%
% kp1 = kp/Np1
% kp2 = kp/Np2
%
% where kp is the stiffness of an individual PEVK residue
%
% Now coming back to our two springs
% l = f1/(kig/Ni1)) + f1/(kp/Np1)
% l = f1*Ni1/kig + f1*Np1/kp
% l = f1*(Ni1*kp + Np1*kig)/(kig*kp)
% l*kig = f1*(Ni1 + Np1*R)
%
% where R = kig/kpevk
%
% and
% l = f2*Ni2/kig + f2*Np2/kp
% l =     f2*(Ni2*kp + Np2*kig)/(kig*kp)
% l*kig = f2*(Ni2 + Np2*R)
%
% Now we're going to massage the previous equation so that we can
% express f2 in terms of f1
%
% l*kig = [f2*(Ni2 + Np2*R)/(Ni1 + Np1*R)]*(Ni1 + Np1*R)
%% 
% Thus
%
% f2*(Ni2 + Np2*R)/(Ni1 + Np1*R) = f1
%
% Now we want to evaluate lp2/lp1. We have
%
% lp1 = (f1/kp1)
%     = f2*(Ni2 + Np2*R)/(Ni1 + Np1*R) / (kp/Np1)
%     = f2*Np1*((Ni2 + Np2*R)/(Ni1 + Np1*R)) / kp
%
% lp2 = (f2/kp2)
%     = f2 / (kp/Np2)
%     = f2*Np2 / (kp)
%
% lp2 / lp1 = (f2* Np2 / (kp)) 
%            /(f2*Np1*((Ni2 + Np2*R)/(Ni1 + Np1*R)) / kp)           
%
%           = (Np2/Np1)((Ni1 + Np1*R)/(Ni2 + Np2*R))
%
% Here 2 is the target and 1 is the reference. So we have
lpevk_div_lpevkRef = ( (nIgRef+nPevkRef*kig_kpevk)...
                       /(nIg + nPevk*kig_kpevk)...
                     )*(nPevk/nPevkRef);

linePevk    = linePevkRef.*lpevk_div_lpevkRef;

%Here is the previous method
linePevk_   = linePevkRef.*(ratioPevkIg/ratioPevkIgRef);

diffWithPreviousMethod = lpevk_div_lpevkRef - (ratioPevkIg/ratioPevkIgRef);

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
