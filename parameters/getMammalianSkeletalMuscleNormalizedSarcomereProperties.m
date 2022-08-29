function [sarcomereProperties] = ...
          getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
            scaleOptimalFiberLength,...
            animalName,...
            normFiberLengthAtOneNormPassiveForce,...
            normPevkToActinAttachmentPoint,...
            fitCrossBridgeStiffnessDampingToKirch199490Hz)  
%%
% This function constructs the normalized lengths of the main components
% of a sarcomere and titin. The data on the relative lengths of actin,
% myosin, and the z-line come from Rassier et al. (Fig 3), Gordon et al., and
% Higuchi et al.
%
% The normalized titin model comes from the values reported by Trombitas et
% al. on human soleus sarcomeres, and from Prado et al. for rabbit psoas.
% I have been unable to find data on the titin segment lengths for a cat
% soleus nor any frog muscles. For now I assume that cat soleus titin looks
% like human soleus titin when normalized by the optimal sarcomere length.
%
% In principle, this titin model, being of normalized
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
% Leonard. Technically it isn't rigid, but its pretty close.
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
%  Higuchi H, Yanagida T, Goldman YE. Compliance of thin filaments in skinned 
%  fibers of rabbit skeletal muscle. Biophysical journal. 1995 Sep 1;
%  69(3):1000-10.
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


animalId = nan;

%disp('Remove this human');
%animalName='rabbit';

if(strcmp(animalName,'cat')==1)
  animalId = 1;
elseif(strcmp(animalName,'human')==1)
  animalId = 2;
elseif(strcmp(animalName,'frog')==1)
  animalId = 3;
elseif(strcmp(animalName,'rabbit')==1)
  animalId = 4;
else  
  assert(0,'animalName must be cat, human, frog, or rabbit');
end


[ geo, ...
  halfMyosinBareLength, ...
  halfMyosinLength,...
  zLineLength,...
  actinLength] = ...
    calcSarcomereFilamentLengthsFromActiveForceLengthKeyPoints(animalId);

%From Fig 3.
optSarcomereLength           = 2*zLineLength + 2*actinLength + 2*halfMyosinBareLength;
normHalfMyosinBareLength     = halfMyosinBareLength/optSarcomereLength; 
normHalfMyosinLength         = halfMyosinLength/optSarcomereLength;
normZLineLength              = zLineLength/optSarcomereLength;                     
normActinLength              = actinLength/optSarcomereLength;
normSarcomereLengthZeroForce = geo(1,1)/optSarcomereLength;

%%
% Human skeletal sarcomere geometry is needed to map Trombitas's 
% experimental data to other animals for which detailed information about
% their titin is not known.
%%
[ geoHuman, ...
  halfMyosinBareLengthHuman, ...
  halfMyosinLengthHuman,...
  zLineLengthHuman,...
  actinLengthHuman] = ...
    calcSarcomereFilamentLengthsFromActiveForceLengthKeyPoints(2);

optSarcomereLengthHuman           = 2*zLineLengthHuman ...
                                    + 2*actinLengthHuman ...
                                    + 2*halfMyosinBareLengthHuman;

normHalfMyosinBareLengthHuman     = halfMyosinBareLengthHuman/optSarcomereLengthHuman; 
normHalfMyosinLengthHuman         = halfMyosinLengthHuman/optSarcomereLengthHuman;
normZLineLengthHuman              = zLineLengthHuman/optSarcomereLengthHuman;                     
normActinLengthHuman              = actinLengthHuman/optSarcomereLengthHuman;
normSarcomereLengthZeroForceHuman = geoHuman(1,1)/optSarcomereLengthHuman;    

                                           
%Check
tol = sqrt(eps);
assert( abs(geo(1,1)-optSarcomereLength*normSarcomereLengthZeroForce) < tol);
assert( abs(geo(1,5)-optSarcomereLength*2*(...
          normZLineLength + normActinLength + normHalfMyosinLength)) < tol);

% The geometry used for the different parts of titin come from Trombitas et
% al's measurements on the titin in a human soleus sample which is known to
% have a long PEVK segment. I haven't yet come across equivalent geometric 
% data on the titin in a cat soleus. I've directly asked W.Herzog, M.DuVall,
% and K.Nishikawa and none of them have seen data on cat soleus titin. So for
% now I'm using the data of Trombitas et al. on the relative elongation of the 
% proximal Ig and PEVK segments for fitting purposes.
%
% Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
% muscle titin revealed by immunolabeling with the anti-titin antibody 9D10. 
% Journal of structural biology. 1998 Jan 1;122(1-2):188-96.


%Reference lengths from human soleus titin (Trombitas et al.)
numDomainsIgPHuman  = 68;
numResiduesPevkHuman= 2174;
numDomainsIgDHuman  = 22;

maxIgDomainStrain_um    = 25/1000;
maxPevkResidueStrain_um = 0.38/1000;
%%
% Titin refrence model: 
%
%   Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
%   muscle titin revealed by immunolabeling with the anti-titin antibody 9D10. 
%   Journal of structural biology. 1998 Jan 1;122(1-2):188-96.
%%

zLineToT12LengthHuman = 0.1; %From Trombitas et al.
zLineToT12Length      = 0.1; %I'm assuming this is the same

flag_assumeUniformIgDomainLength = 1;

fileTrombitas1998Figure5 =...
  ['experiments/TrombitasGreaserFrenchGranzier1998/Trombitas1998_Figure5.csv'];

[lineHumanZToPevkP, lineHumanZToPevkD] = ...
    fitLinearModelToTitinSegmentElongation(fileTrombitas1998Figure5,...
      numDomainsIgPHuman,numResiduesPevkHuman,numDomainsIgDHuman,...
      halfMyosinLengthHuman,zLineToT12LengthHuman,...
      flag_assumeUniformIgDomainLength);

%%
% If the geometry of the animal's titin filament is known, scale the
% elongation data from Trombitas et al. such that 
% 1. The total length of titin is still a half sarcomere
% 2. The elongation of each segment is inversely proportional to its length
%
%%

lineZToPevkP = zeros(2,1);
lineZToPevkD = zeros(2,1);

lTitinHumanOpt =  ...
calcTitinSegmentLengths(  optSarcomereLengthHuman, ...
                          lineHumanZToPevkP, ...
                          lineHumanZToPevkD, ...
                          halfMyosinLengthHuman, ...
                          optSarcomereLengthHuman, ...
                          zLineToT12LengthHuman);

sarcomereLengthHumanOneFpeN = normFiberLengthAtOneNormPassiveForce ...
                                *optSarcomereLengthHuman;

lTitinHumanFisoN =  ...
calcTitinSegmentLengths(  sarcomereLengthHumanOneFpeN, ...
                          lineHumanZToPevkP, ...
                          lineHumanZToPevkD, ...
                          halfMyosinLengthHuman, ...
                          optSarcomereLengthHuman, ...
                          zLineToT12LengthHuman);
                       
fprintf('%e\t%e\tT12\n'          , lTitinHumanOpt.lT12,      lTitinHumanOpt.lT12Norm);
fprintf('%e\t%e\tIgP\n'          , lTitinHumanOpt.lIgp,      lTitinHumanOpt.lIgpNorm);
fprintf('%e\t%e\tPEVK\n'         , lTitinHumanOpt.lPevk,     lTitinHumanOpt.lPevkNorm);
fprintf('%e\t%e\tIgD\n'          , lTitinHumanOpt.lIgdTotal, lTitinHumanOpt.lIgdTotalNorm);
fprintf('%e\t%e\tIgDFixedHuman\n', lTitinHumanOpt.lIgdFixed, lTitinHumanOpt.lIgdFixedNorm);
fprintf('%e\t%e\tIgDFreeHuman\n' , lTitinHumanOpt.lIgdFree,  lTitinHumanOpt.lIgdFreeNorm);

numDomainsIgP     = nan;    
numResiduesPevk   = nan;     
numDomainsIgD     = nan; 

switch animalId
  case 1
    %Cat
    %
    %Scale the data we have from the human soleus by the ratio of
    %optimal sarcomere lengths

    %Since we are assuming that the titin is a proportionate scaling of
    %a human soleus titin, we can calculate the number of prox Ig domains,
    %PEVK residues, and distal Ig domains in our model titin. This is useful
    %later when fitting to data: if the geometry of the prox. Ig and 
    %PEVK residues are unknown. Because this reference is being used for 
    % fitting, I will not round the results even though a fraction of an 
    % Ig domain (or PEVK residue) does not make physical sense.

    numDomainsIgP   =   68*(optSarcomereLength/optSarcomereLengthHuman);    
    numResiduesPevk = 2174*(optSarcomereLength/optSarcomereLengthHuman);     
    numDomainsIgD   =   22*(optSarcomereLength/optSarcomereLengthHuman);     

  case 2
    %Human (soleus titin).
    %
    % Note: soleus titin has a particularly long compliant PEVK segment. 
    %
    %Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
    %muscle titin revealed by immunolabeling with the anti-titin antibody 9D10.
    %Journal of structural biology. 1998 Jan 1;122(1-2):188-96.
    numDomainsIgP     = numDomainsIgPHuman  ;    
    numResiduesPevk   = numResiduesPevkHuman;     
    numDomainsIgD     = numDomainsIgDHuman  ;   

  case 3
    %Frog
    %Scale the data we have from the human soleus by the ratio of
    %optimal sarcomere lengths

    assert(0,['Error: no frog, nor amphibian titin, has had its segment',...
              ' dimensions measured. As a result, here we will not extrapolate',...
              ' results from a mammal to a frog']);

    numDomainsIgP     = nan;    
    numResiduesPevk   = nan;     
    numDomainsIgD     = nan;   

  case 4
     %Rabbit: psoas
     %
     numDomainsIgP   = 50;   
     numResiduesPevk = 800;    
     numDomainsIgD   = 22;
     % For the psoas muscle from a rabbit (Prado et al.). Note that there 
     % are isoforms of titin within a rabbit that have a similar molecular
     % weight as human titin. 
     %
     %Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. 
     % Isoform diversity of giant proteins in relation to passive and active 
     % contractile properties of rabbit skeletal muscles. The Journal of 
     % general physiology. 2005 Nov;126(5):461-80.    
  otherwise
    assert(0,['Error: animalId ',num2str(animalId),' not recognized']);
end


[lineZToPevkP, lineZToPevkD] = ...
    scaleTitinElongationFunction(...
        optSarcomereLength, lTitinHumanOpt.lT12, halfMyosinLength,...
            numDomainsIgP, numResiduesPevk, numDomainsIgD, ...
        optSarcomereLengthHuman, lTitinHumanOpt.lT12, halfMyosinLengthHuman,......
            numDomainsIgPHuman, numResiduesPevkHuman, numDomainsIgDHuman,...
        maxIgDomainStrain_um, maxPevkResidueStrain_um,...
        lineHumanZToPevkP, lineHumanZToPevkD);



dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                      'Sarcomere Length','PEVK Width (um)',...
                      {'a','b'},'Trombitas 1998 Figure 5');  

flag_debugTitinFit =0;
if(flag_debugTitinFit==1)
  figTitinFit =figure;
  for z=1:1:length(dataTrombitas1998Figure5)
    plot(dataTrombitas1998Figure5(z).x,...
         dataTrombitas1998Figure5(z).y,'o','Color',[1,1,1].*0.5,...
         'MarkerFaceColor',[1,1,1].*0.5);
    hold on;
  end
  ls = [2.25;4.75];
  ligp = [ls, ones(size(ls))]*lineHumanZToPevkP;
  plot(ls, ligp,'--','Color',[1,0,0],'DisplayName','Human-Pevk-P');
  hold on;
  lpevk = [ls, ones(size(ls))]*lineHumanZToPevkD;
  plot(ls, lpevk,'--','Color',[1,0,1],'DisplayName','Human-Pevk-P');
  hold on;

  ligp = [ls, ones(size(ls))]*lineZToPevkP;
  plot(ls, ligp,'-','Color',[1,0,0],'LineWidth',1,'DisplayName',[animalName,'-Pevk-P']);
  hold on;
  lpevk = [ls, ones(size(ls))]*lineZToPevkD;
  plot(ls, lpevk,'-','Color',[1,0,1],'LineWidth',1,'DisplayName',[animalName,'-Pevk-P']);
  hold on;  
  xlabel('Sarcomere Length');
  ylabel('Distance from Z-line');
  here=1;
  box off;
  legend('Location','NorthWest');
  legend boxoff;
  
end

%%
% Evaluate the lengths of the segments of titin at optimal, and at
% the normalized length where 1 passive force is reached.
%%

lTitinOpt =  calcTitinSegmentLengths(  ...
                optSarcomereLength, ...
                lineZToPevkP, ...
                lineZToPevkD, ...
                halfMyosinLength, ...
                optSarcomereLength, ...
                zLineToT12Length);


sarcomereLengthOneFpeN = normFiberLengthAtOneNormPassiveForce ...
                        *optSarcomereLength;

lTitinFisoPassive =  calcTitinSegmentLengths(  ...
                  sarcomereLengthOneFpeN, ...
                  lineZToPevkP, ...
                  lineZToPevkD, ...
                  halfMyosinLength, ...
                  optSarcomereLength, ...
                  zLineToT12Length);



stretchHalfLce = 0.5*(max(dataTrombitas1998Figure5(1).x)...
                     -min(dataTrombitas1998Figure5(1).x));

normStretchRateHumanIgP     = lineHumanZToPevkP(1,1);
normStretchRateHumanPevk    = (lineHumanZToPevkD(1,1)-lineHumanZToPevkP(1,1));
normStretchRateHumanIgDFree = (0.5 -(normStretchRateHumanIgP+normStretchRateHumanPevk));

normStretchRateIgP     = lineZToPevkP(1,1);
normStretchRatePevk    = (lineZToPevkD(1,1)-lineZToPevkP(1,1));
normStretchRateIgDFree = (0.5 -(normStretchRateIgP+normStretchRatePevk));
normStretchHalfLce     = stretchHalfLce/optSarcomereLengthHuman;

fprintf('%e\t%e\tNorm. IgP Stretch Rate human vs %s \n'   , ...
    normStretchRateHumanIgP, normStretchRateIgP, animalName);
fprintf('%e\t%e\tNorm. Pevk Stretch Rate human vs %s \n'  , ...
    normStretchRateHumanPevk, normStretchRatePevk, animalName);
fprintf('%e\t%e\tNorm. IgD Stretch Rate human vs %s \n'   , ...
    normStretchRateHumanIgDFree, normStretchRateIgDFree, animalName);
fprintf('%e\tlce stretch (Trobitas et al.) \n'         , ...
    normStretchHalfLce);


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
% I'm setting the titin fraction to be the average of the values
% values reported by Prado which is 44%, with the remaining 56% being due 
% to the ecm.
%
% Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
% diversity of giant proteins in relation to passive and active contractile 
% properties of rabbit skeletal muscles. The Journal of general physiology. 
% 2005 Nov;126(5):461-80.
%

ecmPradoAverage = 0.56;
%Lacking any better data, I'll assign the
%ecmCellularMatrixPassiveForceFraction to be the average across 5 rabbit
%skeletal muscles as reported by Prado et al.

extraCellularMatrixPassiveForceFraction = nan; %1-0.5*(0.24+0.57);

switch animalId
    case 1
        %cat-soleus
        extraCellularMatrixPassiveForceFraction = ecmPradoAverage; 
    case 2
        %human
        extraCellularMatrixPassiveForceFraction = ecmPradoAverage;         
    case 3
        %frog
        assert(0, 'Error: no ECM data available on frog skeletal muscle');
    case 4
        %rabbit
        titinFraction = 0.728-0.158;
        extraCellularMatrixPassiveForceFraction = 1-titinFraction;
        % extraCellularMatrixPassiveForceFraction should probably be 
        % renamed to extraMyofibrillarPassiveForceFraction.
        
    otherwise
        assert(0,'Error: animalId not recognized');
end


    

switch animalId
  case 1

    %cat
    %I'm going to assume that feline titin geometry in skeletal muscle
    %is a scaled version of that of a human soleus. Probably this isn't
    %correct. The alternative is to use the geometry of a  rabbit psoas 
    %titin isoforms measured by Prado et al. This was one of the lightest
    %isoforms measured by Prado et al. (thus having shorter prox. Ig and 
    %PEVK segments).

    lContourIGPNorm     = (68*(25/1000))    / optSarcomereLengthHuman;
    lContourPEVKNorm    = (2174*(0.38/1000))/ optSarcomereLengthHuman;
    lContourIGDFreeNorm = (28*(25/1000))    / optSarcomereLengthHuman; 

  case 2

    % Define the contour lengths of the prox Ig, PEVK, and distal Ig segments   
    %
    % For a human soleus muscle titin's geometry (at least 1 isoform) has been
    % measured by Trombitas
    %
    %   68 prox. Ig domains that can maximally extend to 25 nm (DuVall et al.)
    %   2174 PEVK residues that have a maximum length of 0.38 nm (Cantor & Schimmel)
    %   28 distal Ig domains that can maximally extend to 25 nm 
    %
    %Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
    %muscle titin revealed by immunolabeling with the anti-titin antibody 9D10.
    %Journal of structural biology. 1998 Jan 1;122(1-2):188-96.
    %
    % Cantor CR, Schimmel PR. Biophysical Chemistry, Part I: The Conformation of 
    % Biological Molecules. Journal of Solid-Phase Biochemistry. 1980;5(3).
    %
    % *Note: the 0.38nm is mentioned on 254 as the sum of the bond lengths. In reality
    %        the bond lengths likely cannot be stretched into a line before the
    %        titin filament fails
    % 
    % DuVall MM, Gifford JL, Amrein M, Herzog W. Altered mechanical properties of 
    % titin immunoglobulin domain 27 in the presence of calcium. European Biophysics 
    % Journal. 2013 Apr;42(4):301-7.
    %
    lContourIGPNorm     = (68*(25/1000))    / optSarcomereLengthHuman;
    lContourPEVKNorm    = (2174*(0.38/1000))/ optSarcomereLengthHuman;
    lContourIGDFreeNorm = (28*(25/1000))    / optSarcomereLengthHuman;  

  case 3    

    %Frog
    assert(0,['Error: lacking any estimate of the contour lengths of ',...
                  'the titin segments from frog skeletal muscle titin']);
  case 4


    %%
    % A rabbit psoas has a titin molecule with
    %
    %   50 prox. Ig domains
    %   800 PEVK residues 
    %   22 distal Ig domains 
    %
    % Prado makes it clear that there the range of molecular weights of titin
    % vary quite a bit within the rabbit muscles that were analyzed. Some 
    % of the muscles in a rabbit approach titin in the 3.7 kD range, which 
    % would be consistent with human soleus titin. The size of the PEVK segment
    % seems to be affected most by the molecular weight.
    %
    % This information is really only needed (for the current paper) to 
    % replicate Leonard, Joumaa and Herzog 2010 ... which was performed on a 
    % rabbit psoas muscle.
    %
    %Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
    % diversity of giant proteins in relation to passive and active contractile 
    % properties of rabbit skeletal muscles. The Journal of general physiology. 
    % 2005 Nov;126(5):461-80.
    %
    %%
    
    lContourIGPNorm     = (50*(25/1000))    / optSarcomereLength;
    lContourPEVKNorm    = (800*(0.38/1000))/ optSarcomereLength;
    lContourIGDFreeNorm = (22*(25/1000))    / optSarcomereLength;        

end



% 2022/06/11
% M.Millard
%
titinModelStickySpring      =0; 
%As in the PEVK element viscously sticks to actin
% The idea is similar to that originally proposed by Rode et al., though the 
% details differ: 
%   - here there is one point in the PEVK segment that bonds to actin
%   - the bond is viscous in nature and slides due to imbalanced forces
%
% Rode C, Siebert T, Blickhan R. Titin-induced force enhancement and force 
% depression: a ‘sticky-spring’mechanism in muscle contractions?. Journal of 
% theoretical biology. 2009 Jul 21;259(2):350-60.
%
% Though I've used this model for nearly 2 years, there are two problems with
% it: 
%
% 1. Although I can replicate Leonard, Joumaa, and Herzog using this model, the
%    PEVK + IgD segments are required to flex far far beyond their maximum 
%    contour lengths. 
%
%  2022/06/15
%    W.Herzog mentioned that there is up to 1 um of sarcomere hetereogeneity.
%    This is on the order of magnitude that the PEVK segment could still be
%    within actin overlap.
%
%    I find the sarcomere hetereogeneity argument a bit hollow: to have an
%    average value beyond the total contour length of titin some titins
%    would be shorter (ok) and some would be longer (which means they are
%    broken). After talking to Sven about the practicalities of these
%    experiments he thinks its also possible that the number of sarcomeres
%    were lost track of and/or the sarcomeres on the boundaries were
%    perhaps doing somethin odd.
%
%    Also, part of the odd lengths I observed in my simulations happened because
%    I did not include a sharp increase in force as the Ig and PEVK elements
%    reached their respective contour lengths. I'm updating this now.
%
% 2. This model would not be able to replicate the force
%    traces in Figure 1 of Hisey et al.: when starting from a short length the 
%    PEVK-actin bond would be much closer to the Z-line than the trials that
%    begin at a longer length. Where the biological muscle ends up with a 
%    peak and RFE that is very similar between trials, the simulation would
%    produce very different force profiles: the trial in which the CE was
%    activated at the shortest length would produce far larger forces than the
%    other trials.
%
%  2022/06/16
%    It's possible that if the proximal Ig segment has a finite slack length
%    that the titin-actin attachment point does not change much at short
%    lengths. This is somewhat dependent on what happens to titin as the CE
%    becomes short:
%
%    :If prox. Ig has a slack length and the attachment point is the N2A epitope
%     then the point of attachment will have a lower bound of the slack length
%     of the Ig segment. This is a bit problematic (in the current simulation)
%     because the PEVK segment is about half as stiff as it would need to be
%     to develop the necessary additional tension to replicate Herzog & Leonard
%     2002 with an attachment point at the N2A epitope. It's true that the PEVK 
%     element becomes stiffer with activation, so maybe.
%
%    :If prox and distal Ig segments have a slack length and the attachment 
%     point is in the middle of the PEVK segment, then what happens at short 
%     lengths is not straight forward:
%
%       -If the prox and distal Ig segments remain at their slack lengths and 
%        the PEVK stretches, then at shorter CE lengths the attachment point 
%        will occur at shorter lengths: the model will fail to replicate 
%        Hisey et al. under these circumstances.
%
%    To capture what happens at short lengths I would probably have to introduce
%    another state so that I can track the IgP/PEVK and IgD/PEVK boundaries, as
%    every point in the PEVK segment is a potential binding site to actin. 
%    I think this will best be left for a later model.
%  
% Hisey B, Leonard TR, Herzog W. Does residual force enhancement increase with 
% increasing stretch magnitudes?. Journal of biomechanics. 
% 2009 Jul 22;42(10):1488-92.
%

titinModelActiveSpring        =1; 
% Deprecated: do not use.
%
% After an email exhange with W.Herzog, he mentioned that his hypothesis is that
% the PEVK segment, and perhaps also the Ig segment stiffen under exposure to 
% calcium and this is what produces the large response in force. Single molecule
% studies of both PEVK (Labeit et al.) and the I27 domains of Ig (DuVall) reveal
% that both of these segments increase in stiffness when exposed to calcium. 
% Calcium will bind to the PEVK region specifically (Tatsumi et al.). 
% 
% This model is an extreme version of this phenomenon: here the PEVK
% segment would effectively become rigid when activated. During active
% lengthening the prox/distal Ig segments would lengthen. This mechanism
% would be able to reproduce Hisey et al. However, it is a very strange
% idea. As of 2022/07 I've abandoned it in favour of the sticky-spring
% titin-actin interaction: there is more evidence that the PEVK segment
% bonds, somehow, to actin than there is for the PEVK segment becoming
% rigid. 
%
% Although mechanically realistic models of titin exist (Schappacher-Tilp et al.) these
% require hundreds of state variables to represent the domains. Here we will use
% a very simple lumped model to simulate the PEVK region's resistence to 
% lengthening when exposed to calcium: it will be placed in parallel with a big
% damper which is calcium activated, and saturates under low activation (perhaps
% as low as 1%, Fukutani et al.).
%
% To do: check if pCa4 is consistent with a Ca^2+ concentraction in vivo
%
% Labeit D, Watanabe K, Witt C, Fujita H, Wu Y, Lahmers S, Funck T, Labeit S, 
% Granzier H. Calcium-dependent molecular spring elements in the giant protein 
% titin. Proceedings of the National Academy of Sciences. 2003 Nov 11;100(23):13716-21.
%
% Tatsumi R, Maeda K, Hattori A, Takahashi K. Calcium binding to an elastic 
% portion of connectin/titin filaments. Journal of Muscle Research & Cell 
% Motility. 2001 Feb;22(2):149-62.
%
% DuVall MM, Gifford JL, Amrein M, Herzog W. Altered mechanical properties of titin 
% immunoglobulin domain 27 in the presence of calcium. European Biophysics 
% Journal. 2013 Apr;42(4):301-7.
%
% Schappacher-Tilp G, Leonard T, Desch G, Herzog W. A novel three-filament model 
% of force generation in eccentric contraction of skeletal muscles. PloS one. 
% 2015 Mar 27;10(3):e0117634.
%
% Fukutani A, Herzog W. Residual Force Enhancement Is Preserved for Conditions 
% of Reduced Contractile Force. Medicine and Science in Sports and Exercise. 
% 2018 Jun 1;50(6):1186-91.
%

titinModel = titinModelStickySpring;
%normPevkToActinAttachmentPoint = 0.;
normContourLengthTitinProximal = 0;
normContourLengthTitinDistal = 0;
normLengthTitinFixed = 0;

if(titinModel == titinModelStickySpring)
    normContourLengthTitinProximal = lContourIGPNorm ...
        + normPevkToActinAttachmentPoint*lContourPEVKNorm;
    normContourLengthTitinDistal = (1-normPevkToActinAttachmentPoint)*lContourPEVKNorm ...
        + lContourIGDFreeNorm; 
    normLengthTitinFixed = lTitinOpt.lT12Norm + lTitinOpt.lIgdFixedNorm;
end

if(titinModel == titinModelActiveSpring)
    normContourLengthTitinProximal = lContourIGPNorm; ...
    normContourLengthTitinDistal = lContourPEVKNorm + lContourIGDFreeNorm; 
    normLengthTitinFixed = lTitinOpt.lT12Norm + lTitinOpt.lIgdFixedNorm;
end


%            'IGDFreeNormLengthAtOptimalFiberLengthHuman'    , lTitinHumanOpt.lIgdFreeNorm, ...
%            'IGDFixedNormLengthAtOptimalFiberLengthHuman'   , lTitinHumanOpt.lIgdFixedNorm, ...

disp('  Note: set normECMDamping to 0, from 1e-4');
sarcomereProperties = ...
    struct( 'normActinLength'                   , normActinLength,...  
            'smoothStepFunctionRadius' ,0.01,...
            'normMyosinHalfLength'                          , normHalfMyosinLength,...
            'normMyosinBareHalfLength'                      , normHalfMyosinBareLength,...
            'normZLineLength'                               , normZLineLength,...   
            'normSarcomereLengthZeroForce'                  , normSarcomereLengthZeroForce,...  
            'normFiberLengthAtOneNormPassiveForce'          , normFiberLengthAtOneNormPassiveForce,...
            'ZLineToT12NormLengthAtOptimalFiberLength'      , lTitinOpt.lT12Norm,...
            'IGPNormLengthAtOptimalFiberLength'             , lTitinOpt.lIgpNorm,...              
            'PEVKNormLengthAtOptimalFiberLength'            , lTitinOpt.lPevkNorm,...
            'IGDFreeNormLengthAtOptimalFiberLength'         , lTitinOpt.lIgdFreeNorm, ...
            'IGDFixedNormLengthAtOptimalFiberLength'        , lTitinOpt.lIgdFixedNorm, ...            
            'IGDTotalNormLengthAtOptimalFiberLength'        , lTitinOpt.lIgdTotalNorm, ...    
            'normContourLengthTitinProximal'                , normContourLengthTitinProximal,...
            'normContourLengthTitinDistal'                  , normContourLengthTitinDistal,...
            'normLengthTitinFixed'                          , normLengthTitinFixed,...
            'IGPContourLengthNorm'                          , lContourIGPNorm,...
            'PEVKContourLengthNorm'                         , lContourPEVKNorm,...
            'IGDFreeContourLengthNorm'                      , lContourIGDFreeNorm,...            
            'numberOfIGPDomains'                            , numDomainsIgP,...
            'numberOfPEVKResidues'                          , numResiduesPevk,...
            'numberOfIGDDomains'                            , numDomainsIgD, ...
            'IGPNormStretchRate'                            , normStretchRateIgP,...
            'PEVKNormStretchRate'                           , normStretchRatePevk,...
            'IGDFreeNormStretchRate'                        , normStretchRateIgDFree,...            
            'extraCellularMatrixPassiveForceFraction' , ...
                extraCellularMatrixPassiveForceFraction, ...
            'normCrossBridgeStiffness'         , 48.2636,... 
            'normCrossBridgeDamping'           , 0.1629,...
            'normCrossBridgeCyclingDamping'    , 1.,...
            'normMaxActiveTitinToActinDamping' , 65,... 
            'normPassiveTitinToActinDamping'   , 0.25, ...
            'normPevkToActinAttachmentPoint'   , normPevkToActinAttachmentPoint,... 
            'slidingTimeConstant', 0.001,...
            'forceVelocityCalibrationFactor',1.15,...          
            'activationTimeConstant',0.03,...
            'deactivationTimeConstant',0.03,...         
            'normECMDamping' , 0., ...
            'scaleTitinDistal', 1.0,...
            'scaleECM',  1.0,...
            'scaleTitinProximal',  1.0,...
            'fitCrossBridgeStiffnessDampingToKirch199490Hz',...
             fitCrossBridgeStiffnessDampingToKirch199490Hz,...
             'lowActivationThreshold',0.05,...
             'lowActivationGain',1000,...
             'titinModelType',titinModelStickySpring,...
             'normPassivePevkDamping', 0.25,...
             'normActivePevkDamping',100,...
             'activationThresholdTitin',0.01); 

