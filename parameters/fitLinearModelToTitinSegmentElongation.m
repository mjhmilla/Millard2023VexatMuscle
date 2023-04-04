function [cZToPEVKp, cZToPEVKd] = fitLinearModelToTitinSegmentElongation(...
	fileTrombitas1998Figure5, numDomainsIgP, numResiduesPevk, numDomainsIgD,...
  halfMyosinLength, lT12, flag_assumeUniformIgDomainLength)
%%
% This function fits a line to the digitized data from Trombitas et al. 1998
% Fig.5 and returns the coefficients of the line discribing the distance 
% between Z-line and the proximal end of the PEVK segment, and the distance
% between the Z-line and the distal end of the PEVK segment. To use these
% coefficients to evaluate the length of titin's segments given a sarcomere
% length:
%
%  optSarcomereLengthHuman = 2.73; 
%  lT12          = 0.1; %length between the Z-line and the T12 eptope is 0.1 um
%
%  loptIGP       = [optSarcomereLengthHuman, 1]*cZToPEVKp - lT12;
%  loptPEVK      = [optSarcomereLengthHuman, 1]*cZToPEVKd ...
%                - [optSarcomereLengthHuman, 1]*cZToPEVKp;
%  loptIGDTotal  = 0.5*optSarcomereLengthHuman -(lT12+loptIGP+loptPEVK);
%  loptIGDFixed  = halfMyosinLengthHuman;
%  loptIGDFree   = loptIGDTotal-loptIGDFixed;
% 
%
% Trombitas K, Greaser M, French G, Granzier H. PEVK extension of human soleus 
% muscle titin revealed by immunolabeling with the anti-titin antibody 9D10. 
% Journal of structural biology. 1998 Jan 1;122(1-2):188-96.
%
% @return cZToPEVKp: a 2x1 vector of the offset (cZToPEVKp(2,1)) and slope 
%   (cZToPEVKp(1,1)) that describes how the distance between the Z-line and
%   the proximal end of the PEVK segment varies as the sarcomere is lengthened.
%   The offset and slope assume that the sarcomere length is given in 
%   micrometers.
%
% @return cZToPEVKd: a 2x1 vector of the offset (cZToPEVKd(2,1)) and slope 
%   (cZToPEVKd(1,1)) that describes how the distance between the Z-line and
%   the distal end of the PEVK segment varies as the sarcomere is lengthened.
%   The offset and slope assume that the sarcomere length is given in 
%   micrometers.
%
%
%%

dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                      'Sarcomere Length','PEVK Width (um)',...
                      {'a','b'},'Trombitas 1998 Figure 5');  

%Least squares fit
%Z-to-Igp
ba = dataTrombitas1998Figure5(1).y;
Aa = [dataTrombitas1998Figure5(1).x, ones(size(ba))];
cZToPEVKp = pinv(Aa'*Aa)*(Aa'*ba);

%lT12 = 0.1;
%l

%Z-to-Pevk
bb = dataTrombitas1998Figure5(2).y;
Ab = [dataTrombitas1998Figure5(2).x, ones(size(bb))];
cZToPEVKd = pinv(Ab'*Ab)*(Ab'*bb);

cZToPEVKpAdj = cZToPEVKp;
cZToPEVKdAdj = cZToPEVKd;  

if(flag_assumeUniformIgDomainLength==1)
  linePevk     = cZToPEVKd - cZToPEVKp;
  lineIgTotal = [0.5;0] - linePevk ...
                 -[0;halfMyosinLength] - [0;lT12];
%numDomainsIgP, numResiduesPevk, numDomainsIgD
  lineIgP = lineIgTotal.*(numDomainsIgP/(numDomainsIgP+numDomainsIgD));
  lineIgD = lineIgTotal.*(numDomainsIgD/(numDomainsIgP+numDomainsIgD));

  cZToPEVKpAdj = lineIgP + [0;lT12];
  cZToPEVKdAdj = cZToPEVKpAdj + linePevk;  

end


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
  ligp = [ls, ones(size(ls))]*cZToPEVKp;
  plot(ls, ligp,'--','Color',[1,0,0]);
  hold on;
  lpevk = [ls, ones(size(ls))]*cZToPEVKd;
  plot(ls, lpevk,'--','Color',[1,0,1]);
  hold on;

  ligp = [ls, ones(size(ls))]*cZToPEVKpAdj;
  plot(ls, ligp,'-','Color',[1,0,0],'LineWidth',1);
  hold on;
  lpevk = [ls, ones(size(ls))]*cZToPEVKdAdj;
  plot(ls, lpevk,'-','Color',[1,0,1],'LineWidth',1);
  hold on;  
  xlabel('Sarcomere Length');
  ylabel('Distance from Z-line');
  here=1;
  
end

cZToPEVKp = cZToPEVKpAdj;
cZToPEVKd = cZToPEVKdAdj;  