function [sarcomereProperties] = ...
          getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
            scaleOptimalFiberLength,...
            flag_Cat1_Human2,...
            fitCrossBridgeStiffnessDampingToKirch199490Hz)  
%%
% This function constructs the normalized lengths of the main components
% of a sarcomere and titin. The data on the relative lengths of actin,
% myosin, and the z-line come from Rassier et al. (Fig 3) and Gordon et al.
%
% The normalized titin model comes from the values reported by Trombitas et
% al. on human soleus sarcomeres. This titin model, being of normalized
% length can be applied to any of the sarcomere models. It is assumed that
% the distal Ig segment that overlaps myosin at the optimal fiber length is
% permanently bound to myosin and thus is not free to flex.
%
% The fraction of the passive force length curve that is attributed to
% the extracellular matrix comes from the average of the maximum and
% minimum fraction observed by Prado et al. in the study they conducted
% on rabbit skeletal muscle.
%
% Finally, the terminal stiffness of the section of titin between the
% N2A element and the distal part of Ig region that is bound to myosin is
% estimated using the data from the first trace of Fig. 7A of Herzog and
% Leonard.
%
% References                    
%  Rassier DE, MacIntosh BR, Herzog W. Length dependence of active force 
%  production in skeletal muscle. Journal of applied physiology. 
%  1999 May 1;86(5):1445-57.
%
%  Gordon, A. M., A. F. Huxley, and F. J. Julian. The variation in
%  isometric tension with sarcomere length in vertebrate muscle
%  fibres. J. Physiol. (Lond.) 184: 170–192, 1966. 
%
%  Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
%  muscle titin revealed by immunolabeling with the anti-titin antibody 9D10. 
%  Journal of structural biology. 1998 Jan 1;122(1-2):188-96.
%
%  Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
%  diversity of giant proteins in relation to passive and active contractile 
%  properties of rabbit skeletal muscles. The Journal of general physiology. 
%  2005 Nov;126(5):461-80.
%
%  Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
%  muscle: a new mechanism. Journal of Experimental Biology. 2002 May 1;
%  205(9):1275-83.
%
%
%%

geoCat   = [1.27, 1.7, 2.34, 2.51, 3.94];
geoHuman = [1.27, 1.7, 2.64, 2.81, 4.24];

%The frog geometry has been blocked: the titin geometry below comes from
%a human soleus. Perhaps it is similar to cat skeletal muscel. However,
%there is no reason (that I) have to expect that it is similar to the titin
%molecule in frog skeletal muscle.
%
%geoFrog  = [1.27, 1.7,  2.0,  2.2, 3.6 ];
geo = geoCat;
switch flag_Cat1_Human2
  case 1
    geo = geoCat;    
  case 2
    geo = geoHuman;    
  otherwise
    assert(0, 'Error: flagFrog0Cat1Human1 incorrectly set');
end
    
%Note that the length of the bare section is 1/4 of the length
%of the plateau. Why? The first half comes from the fact that 
%the section from geo(1,4) to geo(1,3) is a result of 2 actins moving
%across the flat section. We take half of that again because we're
%interested in half of the bare length.
halfMyosinBareLength =     ( geo(1,4) - geo(1,3) )*0.25;

halfMyosinLength     = 0.5*(geo(1,5)-geo(1,4)) ...
                      +halfMyosinBareLength;                    
zLineLength = 0.05;                    
actinLength          = 0.5*( (geo(1,5))...
                        -2*halfMyosinLength ...
                        -2*zLineLength);
                      
                      

%From Fig 3.
optSarcomereLength   = 2*zLineLength + 2*actinLength;

%0.5*(geo(1,3)+geo(1,4));  

%From Fig 2. - taking the bare length from the frog data reported by Gordon
normHalfMyosinBareLength = halfMyosinBareLength/optSarcomereLength; 

%From Fig 3: length where overlap is lost less optimal
normHalfMyosinLength = halfMyosinLength/optSarcomereLength;

%From Fig 2. - taking the z-line length from the frog data reported by Gordon
normZLineLength = zLineLength/optSarcomereLength; %Rassier et al.                     


normActinLength = actinLength/optSarcomereLength;

%From Fig 3, which is an extrapolation
normSarcomereLengthZeroForce = geo(1,1)/optSarcomereLength;

%Check
tol = sqrt(eps);
assert( abs(geo(1,1)-optSarcomereLength*normSarcomereLengthZeroForce) < tol);
assert( abs(geo(1,5)-optSarcomereLength*2*(...
          normZLineLength + normActinLength + normHalfMyosinLength)) < tol);

% The geometry used for the different parts of titin come from Trombitas et
% al's measurements on the titin in a human soleus sample. I haven't yet
% come across equivalent geometric data on the titin in a cat soleus. 
% This should be a decent approximation because both samples come from
% mammalian skeleltal muscle (the soleus). On the other hand there are many 
% titin isoforms: it is possible there are big differences. There's nothing 
% I can do at this moment to improve these geometric estimates.
%
% Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
% muscle titin revealed by immunolabeling with the anti-titin antibody 9D10. 
% Journal of structural biology. 1998 Jan 1;122(1-2):188-96.

fileTrombitas1998Figure5 = 'experiments/TrombitasGreaserFrenchGranzier1998/Trombitas1998_Figure5.csv';
dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                      'Sarcomere Length','PEVK Width (um)',...
                      {'a','b'},'Trombitas 1998 Figure 5');  

%Least squares fit
%Z-to-Igp
ba = dataTrombitas1998Figure5(1).y;
Aa = [dataTrombitas1998Figure5(1).x, ones(size(ba))];
xa = (Aa'*Aa)\(Aa'*ba);

%lT12 = 0.1;
%l

%Z-to-Pevk
bb = dataTrombitas1998Figure5(2).y;
Ab = [dataTrombitas1998Figure5(2).x, ones(size(bb))];
xb = (Ab'*Ab)\(Ab'*bb);

flag_debugTitinFit =1;
if(flag_debugTitinFit==1)
  figTitinFit =figure;
  for z=1:1:length(dataTrombitas1998Figure5)
    plot(dataTrombitas1998Figure5(z).x,...
         dataTrombitas1998Figure5(z).y,'o','Color',[1,1,1].*0.5,...
         'MarkerFaceColor',[1,1,1].*0.5);
    hold on;
  end
  ls = [2.25;4.75];
  ligp = [ls, ones(size(ls))]*xa;
  plot(ls, ligp,'Color',[1,0,0]);
  hold on;
  lpevk = [ls, ones(size(ls))]*xb;
  plot(ls, lpevk,'Color',[1,0,1]);
  hold on;
  xlabel('Sarcomere Length');
  ylabel('Distance from Z-line');
  here=1;
  
end
                    
%Note All IGP, IGD, PEVK, and T12 lengths come from Trombitas et al. 1997
%     by assuming lopt is 2.725 um, and that the passive element develops
%     1 maximum-isometric-force passively at 1.6 lopt.
optSarcomereLengthHuman = 0.5*(geoHuman(1,3)+geoHuman(1,4));  

halfMyosinBareLengthHuman =( geoHuman(1,4) - geoHuman(1,3) )*0.25;
halfMyosinLengthHuman     = 0.5*(geoHuman(1,5)-geoHuman(1,4)) ...
                           +halfMyosinBareLengthHuman;  

normHalfMyosinLengthHuman     = halfMyosinLengthHuman/optSarcomereLengthHuman;

lT12     = 0.1;
loptIGP  = [optSarcomereLengthHuman, 1]*xa - lT12;
loptPEVK = [optSarcomereLengthHuman, 1]*xb - [optSarcomereLengthHuman, 1]*xa;
loptIGDTotal = 0.5*optSarcomereLengthHuman -(lT12+loptIGP+loptPEVK);
loptIGDFixed = halfMyosinLengthHuman;
loptIGDFree  = loptIGDTotal-loptIGDFixed;



loptT12Norm  = lT12/optSarcomereLengthHuman;
loptIGPNorm  = loptIGP/optSarcomereLengthHuman; %Should be renamed IGP: proximal Ig
loptPEVKNorm = loptPEVK/optSarcomereLengthHuman;
loptIGDTotalNorm  = 1*0.5-(loptT12Norm+loptIGPNorm+loptPEVKNorm); %0.5: applies to a half sarcomere






%IGD: This is the small section of the distal Ig that section that
% is still free to stretch. The rest is bound to myosin and does not 
% stretch at all.

disp('Reminder: remove all lfiso variables');
lfisoIGPNorm  = 0.45812/optSarcomereLengthHuman;
lfisoPEVKNorm = 0.6454/optSarcomereLengthHuman;
lfisoIGDTotalNorm  = (geoHuman(1,5)*0.5/optSarcomereLengthHuman ...
               -(loptT12Norm + lfisoIGPNorm+lfisoPEVKNorm));
                       
%Now we have a normalized titin model. Since the relative length of myosin
%in the cat, human, and frog are different we must update the length of
%IGD. First we set these lengths for a human title model (soleus from
%Trombitas et al) as a reference, and next we set the lengths for the
%animal model.
loptIGDFreeNormHuman   = loptIGDTotalNorm - normHalfMyosinLengthHuman; 
loptIGDFixedNormHuman  = normHalfMyosinLengthHuman;

fprintf('%e\t%e\tT12\n'     , lT12, loptT12Norm);
fprintf('%e\t%e\tIgP\n'     , loptIGP, loptIGPNorm);
fprintf('%e\t%e\tPEVK\n'    , loptPEVK, loptPEVKNorm);
fprintf('%e\t%e\tIgD\n'     , loptIGDTotal, loptIGDTotalNorm);
fprintf('%e\t%e\tIgDFixedHuman\n', loptIGDFixed, loptIGDFixedNormHuman);
fprintf('%e\t%e\tIgDFreeHuman\n' , loptIGDFree, loptIGDFreeNormHuman);

%If the entire distal Ig region overlaps with myosin then the entire
%length of the distal Ig region is bound to myosin.
if(loptIGDFreeNormHuman <= 0)
  loptIGDFreeNormHuman  = 0;
  loptIGDFixedNormHuman = loptIGDTotalNorm;
  disp('Warning: Distal Ig length is completely bound to myosin');  
end

%Repeat for the animal model
loptIGDFreeNorm   = loptIGDTotalNorm - normHalfMyosinLength; 
loptIGDFixedNorm  = normHalfMyosinLength;

%If the entire distal Ig region overlaps with myosin then the entire
%length of the distal Ig region is bound to myosin.
if(loptIGDFreeNorm <= 0)
  loptIGDFreeNorm  = 0;
  loptIGDFixedNorm = loptIGDTotalNorm;
  disp('Warning: Distal Ig length is completely bound to myosin');  
end


lfisoIGDFreeNorm  = lfisoIGDTotalNorm - loptIGDFixedNorm;

stretchHalfLce = 0.5*(max(dataTrombitas1998Figure5(1).x)-min(dataTrombitas1998Figure5(1).x));

normStretchRateIgP     = xa(1,1);
normStretchRatePevk    = (xb(1,1)-xa(1,1));
normStretchRateIgDFree = (0.5 -(normStretchRateIgP+normStretchRatePevk));
normStretchHalfLce     = stretchHalfLce/optSarcomereLengthHuman;

loptIGDFixedCat = normHalfMyosinLength*optSarcomereLength;
loptIGDFreeCat =  loptIGDFreeNorm*optSarcomereLength;


fprintf('%e\t%e\tIgDFixedCat\n'            , loptIGDFixedCat, loptIGDFixedNorm);
fprintf('%e\t%e\tIgDFreeCat\n'             , loptIGDFreeCat, loptIGDFreeNorm);
fprintf('%e\tNorm. IgP Stretch Rate\n'     , normStretchRateIgP);
fprintf('%e\tNorm. Pevk Stretch Rate\n'    , normStretchRatePevk);
fprintf('%e\tNorm. IgD Stretch Rate\n'     , normStretchRateIgDFree);
fprintf('%e\tlce Stretch Rate\n'     , normStretchHalfLce);

here=1;


%
%Page 472 column 2 last paragraph of Prado et al. reports:
%
%     These results show that titin’s relative contribution to
%     total passive stiffness is much higher in some muscles,
%     like psoas (57%) and diaphragm (56%), than in oth-
%     ers, like soleus (24%), EDL (42%), and gastrocne-
%     mius (41%)
%
% Lacking any real information about the titin fraction in a cat's soleus
% I'm setting the titin fraction to be the average of the lowest and highest
% values reported by Prado: (0.24+0.57)*0.5
%
% Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
% diversity of giant proteins in relation to passive and active contractile 
% properties of rabbit skeletal muscles. The Journal of general physiology. 
% 2005 Nov;126(5):461-80.
%
extraCellularMatrixPassiveForceFraction = 0.56; %1-0.5*(0.24+0.57);

% Here I am assuming that the large long positive ramp in the force trace
% of Fig. 7A is entirely due to the PEVK and distal IG region stretching. Thus
% to estimate the stiffness of the PEVK and distal IG region I can get the 
% slope of this line. 
%
% As a note this is intended as a rough initial estimate that later is fitted to
% data.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 2002 May 1;
% 205(9):1275-83.
%                       

fitData =csvread('experiments/HerzogLeonard2002/data/dataHerzogLeonard2002Figure6and7A.csv',1,0);

%Extracted from the musculotendon data on Herzog & Leonard 2002
%See function getFelineSoleusMusculotendonProperties for details.
%Here the values are hard coded because we are using Herzog & Leonard 2002
%to evaluate the normalized stiffness of the PEVK + IGD sections
fiso    = 2.141875213041221e+01;
lopt    = 4.285714285714286e-02 * scaleOptimalFiberLength;
alphaOpt= 1.221730476396031e-01;


%Note: 'F43', 'L43', etc are the column headings in the file.
idxFitF43=2;   % Length 0mm, Activation 0-1-0 
idxFitL43=3;
idxFitF44=4;   % Length 9mm, Activation 0-1-0
idxFitL44=5;
idxFitF45=6;   % Length 6-9mm, Activation 0-1-0
idxFitL45=7;
idxFitF47=8;   % Length 3-9mm, Activation 0-1-0
idxFitL47=9;
idxFitF49=10;  % Length 0-9mm, Activation 0 (Passive)
idxFitL49=11;
idxFitF50=12;  % Length 0-9mm, Activation 0-1-0
idxFitL50=13;

idxOpt= 230;
idx0  = 297; %Beginning of the long active increase in force with length
idx1  = 507; %End " ... "


dlAT0 = fitData(idx0,idxFitL50);
dlAT1 = fitData(idx1,idxFitL50);
fpeAT0 = 0;
fpeAT1 = 0;

idxMax = 1500;

%Get the passive forces at these two lengths
dlerr0=Inf;
dlerr1=Inf;
for i=1:1:idxMax
  if( abs(dlAT0-fitData(idx0,idxFitL49)) < dlerr0 )
    dlerr0 = abs(dlAT0-fitData(idx0,idxFitL49));
    fpeAT0 = fitData(idx0,idxFitF49);
  end
  if( abs(dlAT1-fitData(idx1,idxFitL49)) < dlerr1 )
    dlerr1 = abs(dlAT1-fitData(idx1,idxFitL49));
    fpeAT1 = fitData(idx1,idxFitF49);
  end  
end
%Subtract off the component of the measured force that is due to 
%the extracellular matrix using the estimated fraction of ECM from Prado et al.
extraCellularMatrixPassiveForceFractionHerzogLeonard2002 = 1-0.5*(0.24+0.57);
fAT0 =  fitData(idx0,idxFitF50) ...
       -fpeAT0*extraCellularMatrixPassiveForceFractionHerzogLeonard2002;
fAT1 = fitData(idx1,idxFitF50) ...
      -fpeAT1*extraCellularMatrixPassiveForceFractionHerzogLeonard2002;
mm2m = 0.001;

%Get the fiber length and pennation angle at the first point
l0AT   =  lopt*cos(alphaOpt) + fitData(idx0,idxFitL50).*mm2m; 
fibKin = calcFixedWidthPennatedFiberKinematics(l0AT,0,lopt,alphaOpt);
l0     = fibKin.fiberLength;
alpha0 = fibKin.pennationAngle;

%Get the fiber length and pennation angle at the end point
l1AT   = lopt*cos(alphaOpt) + fitData(idx1,idxFitL50).*mm2m;
fibKin = calcFixedWidthPennatedFiberKinematics(l1AT,0,lopt,alphaOpt);
l1     = fibKin.fiberLength;
alpha1 = fibKin.pennationAngle;

f0  = fAT0/cos(alpha0);
f1  = fAT1/cos(alpha1);

l0N = l0/lopt;
l1N = l1/lopt;
f0N = f0/fiso;
f1N = f1/fiso;

normActiveFiberForce  = [f0N,f1N];
normFiberLength       = [l0N,l1N];

%Now estimate the loss in force due to the active-force-length curve
%Make sure we're on the descending limb 
assert(l0N >= 1. && l1N > l0N);

%From Fig 3.
optSarcomereLengthCat   = geoCat(1,4);%0.5*(geoCat(1,3)+geoCat(1,4));  

%From Fig 3: length where overlap is lost less optimal
normHalfMyosinLengthCat = 0.5*(geoCat(1,5) - optSarcomereLengthCat) ...
                       /optSarcomereLengthCat; 

%From Fig 2. - taking the bare length from the frog data reported by Gordon
normHalfMyosinBareLengthCat = ((geoCat(1,4)-geoCat(1,3))/optSarcomereLengthCat)*0.5; 


%From Fig 2. - taking the z-line length from the frog data reported by Gordon
normZLineLengthCat      = 0.05/optSarcomereLengthCat; %Rassier et al.                     

normActinLengthCat = 0.5 - (normHalfMyosinBareLengthCat+normZLineLengthCat);


falAtlceOpt  = 1;
falAtlceMax  = 0;
lceOpt       = 1;
lceMax       =  2*(normHalfMyosinLengthCat...
                  +normHalfMyosinBareLengthCat...
                  +normActinLengthCat...
                  +normZLineLengthCat);
Dfal_DlceN = (falAtlceOpt-falAtlceMax)/(lceOpt - lceMax);

disp('Remove kNormPevkIgd');
kNormPevkIGd = diff(normActiveFiberForce)/diff(normFiberLength) ...
               - Dfal_DlceN;

% Solve for the stiffness of the Igp region by using the ratio of the
% relative stretches between these regions as observed by Trombitas

deltaIgp     = lfisoIGPNorm      - loptIGPNorm;
deltaIgdFree = lfisoIGDFreeNorm  - loptIGDFreeNorm;    
deltaPEVK    = lfisoPEVKNorm     - loptPEVKNorm;

disp('Remove kNormIgd');
kNormIgp = kNormPevkIGd*(deltaIgdFree+deltaPEVK)/(deltaIgp);
             


disp('  Note: set normECMDamping to 0, from 1e-4');
sarcomereProperties = ...
    struct( 'normActinLength'                   , normActinLength,...  
            'normMyosinHalfLength'              , normHalfMyosinLength,...
            'normMyosinBareHalfLength'          , normHalfMyosinBareLength,...
            'normZLineLength'                   , normZLineLength,...   
            'normSarcomereLengthZeroForce'      , normSarcomereLengthZeroForce,...       
            'ZLineToT12NormLengthAtOptimalFiberLength', loptT12Norm,...
            'IGPNormLengthAtOptimalFiberLength'       , loptIGPNorm,...  
            'PEVKNormLengthAtOptimalFiberLength'      , loptPEVKNorm,...
            'IGDFreeNormLengthAtOptimalFiberLengthHuman'   , loptIGDFreeNormHuman, ...
            'IGDFixedNormLengthAtOptimalFiberLengthHuman'  , loptIGDFixedNormHuman, ...
            'IGDFreeNormLengthAtOptimalFiberLength'        , loptIGDFreeNorm, ...
            'IGDFixedNormLengthAtOptimalFiberLength'       , loptIGDFixedNorm, ...            
            'IGDTotalNormLengthAtOptimalFiberLength'       , loptIGDTotalNorm, ...               
            'PEVKIGDNormStiffness'                    , kNormPevkIGd,...    
            'IGPNormStiffness'                        , kNormIgp,...
            'IGPNormStretchRate'                      , normStretchRateIgP,...
            'PEVKNormStretchRate'                     , normStretchRatePevk,...
            'IGDFreeNormStretchRate'                  , normStretchRateIgDFree,...            
            'extraCellularMatrixPassiveForceFraction' , ...
                extraCellularMatrixPassiveForceFraction, ...
            'normCrossBridgeStiffness'         , 48.2636,... 
            'normCrossBridgeDamping'           , 0.1629,...
            'normCrossBridgeCyclingDamping'    , 1.,...
            'normMaxActiveTitinToActinDamping' , 75,... 
            'normPassiveTitinToActinDamping'   , 0.25, ...
            'normPevkToActinAttachmentPoint'   , 0.6,... 
            'slidingTimeConstant', 0.001,...
            'forceVelocityCalibrationFactor',1.15,...          
            'activationTimeConstant',0.03,...
            'deactivationTimeConstant',0.03,...         
            'normECMDamping' , 0., ...
            'scalePEVK', 1.0,...
            'scaleECM',  1.0,...
            'scaleIGP',  1.0,...
            'fitCrossBridgeStiffnessDampingToKirch199490Hz',...
             fitCrossBridgeStiffnessDampingToKirch199490Hz,...
             'lowActivationThreshold',0.05);
