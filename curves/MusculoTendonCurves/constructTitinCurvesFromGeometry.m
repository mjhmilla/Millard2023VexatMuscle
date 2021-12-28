function [forceLengthIgpCurve, forceLengthIgpInverseCurve,...
          forceLengthPevkCurve, forceLengthPevkInverseCurve] = ...
    constructTitinCurvesFromGeometry(...
                        fiberForceLengthCurve,...                                   
                        forceLengthCurveSettings,...
                        forceLengthECMHalfCurve,...
                        musculotendonProperties,...
                        sarcomereProperties,...
                        muscleName,...
                        flag_computeCurveIntegrals,...
                        flag_useElasticIgD,...
                        flag_useHumanIgDGeometry,...
                        flag_useOctave)

%Cantor CR, Schimmel PR. Biophysical Chemistry, Part I: The Conformation of 
% Biological Molecules. Journal of Solid-Phase Biochemistry. 1980;5(3).
%
%DuVall MM, Gifford JL, Amrein M, Herzog W. Altered mechanical properties 
% of titin immunoglobulin domain 27 in the presence of calcium. European 
% Biophysics Journal. 2013 Apr 1;42(4):301-7.
%
%Trombitás K, Greaser M, Labeit S, Jin JP, Kellermayer M, Helmes M, 
% Granzier H. Titin extensibility in situ: entropic elasticity of 
% permanently folded and permanently unfolded molecular segments. The 
% Journal of cell biology. 1998 Feb 23;140(4):853-9.


forceLengthIgpCurve         =   []; 
forceLengthIgpInverseCurve  =   [];
forceLengthPevkCurve        =   [];
forceLengthPevkInverseCurve =   [];

numberOfProximalIgDomains = 66; %Ig-like domains
numberOfDistalIgDomains   = 28; %Ig/Fn-like domains
%From table 1: foot note of Trombitas et al.

maximumIgDomainLength     = 25.25e-9; 
%Taking the maximum domain length in the presence of Ca
%(24.76 nm ± 0.13 nm) and calcium (25.25 nm ± 0.13 nm)
%From DuVall et al.

contourLengthProximalIg   = numberOfProximalIgDomains*maximumIgDomainLength;
contourLengthDistalIg     = numberOfDistalIgDomains*maximumIgDomainLength;
persistenceLengthIg       = 15e-9; %nm From Trombitás et al.

numberOfPevkResidues      = 2174;
maximumPevkResidueLength  = 0.38e-9; %nm
persistenceLengthPevk     = 2e-9; %nm From Trombitás et al.

%0.32nm
%from Trombitás et al.

%0.38nm
%From Cantor and Schimmel.

contourLengthPevk         = numberOfPevkResidues*maximumPevkResidueLength;
%nm From Trombitás et al.

T = 38+273.15; %Body temperature in Kelvin
kB= 1.38064852e-23;% m2 kg s-2 K-1; %Boltzmann's constant
AIg     = persistenceLengthIg;
APevk   = persistenceLengthPevk;

Tkb_div_A_Ig    = T*(kB/AIg);
Tkb_div_ASq_Ig  = Tkb_div_A_Ig/AIg;

Tkb_div_A_Pevk    = T*(kB/APevk);
Tkb_div_ASq_Pevk  = Tkb_div_A_Pevk/APevk;

sarcomereHexagonSideLength = 44e-9; %Nm.
%From Fig. 2 of
% Huxley AF. Muscle structure and theories of contraction. Prog. Biophys. 
% Biophys. Chem. 1957;7:255-318.

sarcomereHexagonArea = (3*sqrt(3)/2)*(sarcomereHexagonSideLength)^2;
sarcomereMaxActiveStress = 200e3; %200 kPA
sarcomereMaxActiveForce = sarcomereMaxActiveStress*sarcomereHexagonArea;

fmax=(sarcomereMaxActiveForce/6)*5.1;
fmin=0;

% fmax = calcWormLikeChainModelDer(0.95*contourLengthPevk,...
%                                               contourLengthPevk,...
%                                               Tkb_div_A_Pevk,...
%                                               Tkb_div_ASq_Pevk,...
%                                               [0,0,0]);
% 
% fmin = calcWormLikeChainModelDer(0.*contourLengthPevk,...
%                                               contourLengthPevk,...
%                                               Tkb_div_A_Pevk,...
%                                               Tkb_div_ASq_Pevk,...
%                                               [0,0,0]);

npts=1000;
f = [fmin:(fmax-fmin)/(npts-1):fmax]';

lIgP    = zeros(size(f));
lPevk   = zeros(size(f));
lIgD    = zeros(size(f));

for i=2:1:npts
    lIgP(i,1)    = calcWormLikeChainModelInvDer(...
                f(i,1),lIgP(i-1,1),contourLengthProximalIg,Tkb_div_A_Ig,...
                Tkb_div_ASq_Ig,[0,0,0]);
    lIgD(i,1)    = calcWormLikeChainModelInvDer(...
                f(i,1),lIgD(i-1,1),contourLengthDistalIg,Tkb_div_A_Ig,...
                Tkb_div_ASq_Ig,[0,0,0]);
    lPevk(i,1)   = calcWormLikeChainModelInvDer(...
                f(i,1),lPevk(i-1,1),contourLengthPevk,Tkb_div_A_Pevk,...
                Tkb_div_ASq_Pevk,[0,0,0]);
end

optimalSarcomereLength=2.725e-6;

l_Z_T12 = ...
  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength...
  *optimalSarcomereLength;

l_halfMyosin = sarcomereProperties.normMyosinHalfLength ...
    *optimalSarcomereLength;

lce = (l_Z_T12 + lIgP + lPevk + lIgD + l_halfMyosin).*2;
fce = f.*6; %6 titins per half sarcomere.

l_z_IgP     = (l_Z_T12 + lIgP);
l_z_Pevk    = (l_Z_T12 + lIgP + lPevk);




fceNorm = fce ./sarcomereMaxActiveForce;
lceNorm = lce ./ optimalSarcomereLength;

fig = figure;
subplot(1,3,1);
    plot(lceNorm, fceNorm,'b');
    xlabel('Norm. Length (m/m)')
    ylabel('Norm Force (N/N)')
    title('Titin Force-Length Curve (sarcomere-level)')

subplot(1,3,2);
    plot(lIgP./optimalSarcomereLength, ...
         fceNorm,'r','DisplayName','IgP');
    hold on;
    plot(lPevk./optimalSarcomereLength, ...
         fceNorm,'b','DisplayName','Pevk');
    hold on;    
    plot(lIgD./optimalSarcomereLength, ...
         fceNorm,'m','DisplayName','IgD');
    hold on;

    plot((lIgD+lPevk)./optimalSarcomereLength, ...
         fceNorm,'k','DisplayName','IgD');
    hold on;
    
    legend;
    xlabel('Norm. Length (m/m)')
    ylabel('Norm Force (N/N)')
    title('Titin Force-Length Curve')

subplot(1,3,3);
    plot(lce, l_z_IgP,'r','DisplayName','Z-IgP/Pevk');
    hold on;
    plot(lce, l_z_Pevk,'b','DisplayName','Z-Pevk/IgD');
    hold on;
    xlabel('Segment Length (um)');
    ylabel('Sarcomere Length (um)');

    here=1;
