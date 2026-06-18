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

clc;
close all;
clear all;

addpath('parameters/');
addpath('postprocessing/');

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);


% Huxley & Simmons 1981: number crossbridges per half-sarcomere
%2*(700/14.3); %pair of xbridges every 14.3 nm for 700 nm

%Anthony Hessel informed me of some literature that shows that the number
% of crossbridges per half sarcomere is closer to 300! That's a lot
% different than Huxley and Simmons 1981, and of course affects our
% estimates of stiffness.
% 
% Piazzesi G, Reconditi M, Linari M, Lucii L, Bianco P, Brunello E, 
% Decostre V, Stewart A, Gore DB, Irving TC, Irving M. Skeletal muscle 
% performance determined by modulation of number of myosin motors rather 
% than motor force or stroke size. Cell. 2007 Nov 16;131(4):784-95.
%
numCrossbridgesHalfSarcomere = 294*2; %2 cross-bridges per myosin motor


%Linari M, Dobbie I, Reconditi M, Koubassova N, Irving M, Piazzesi G, 
% Lombardi V. The stiffness of skeletal muscle in isometric contraction and 
% rigor: the fraction of myosin heads bound to actin. Biophysical journal. 
% 1998 May 1;74(5):2459-73.
maxActivationProportionOfAttachedCrossbridges = 0.43;

numberOfAttachedCrossbridges=...
  maxActivationProportionOfAttachedCrossbridges*numCrossbridgesHalfSarcomere;

%Units:
% nm: length
% pN: force
projectFolders.output
load(fullfile(projectFolders.output_structs_FittedModels,'tunedRabbitPsoasFibril.mat'));

kTitinProx = tunedRabbitPsoasFibril.curves.forceLengthProximalTitinCurve.dydxEnd(1,2);
kTitinDistal = tunedRabbitPsoasFibril.curves.forceLengthProximalTitinCurve.dydxEnd(1,2);

%Stiffness of the passive titin model
kTitinPassive = 1/((1/kTitinDistal)+(1/kTitinProx));
kTitinActive  = kTitinDistal;

%When activated (in our model) only the distal segment of titin lengthens
ratio_kTitinActive_kTitinPassive = kTitinActive/kTitinPassive;

% Sources
% Viegel et al. 1998: crossbridgeStiffness
% Higuchi et al 1995: actinStiffness, actinLength (rabbit)
% Tajima et al. 1994: myosin stiffness
%  Note this interpretation is that myosin has the same stiffness
%  of 46-68 pN/nm. Since myosin is shorter than actin, this interpretation
%  differs from assuming that both filaments have the same material
%  stiffness.
% Kellermayer et al. 1997: stiffness of titin (rabbit)
sarcomereProperties=struct(...
    'crossbridgeStiffness',[(0.69-0.47),(0.69+0.47)],...
    'actinStiffness',[46,68],...
    'myosinStiffness',[46,68],...
    'myosinLength', 800,...
    'actinLength', 1120);

lengthMLineToMyosinMeanAttachmentPoint = 450; %nm;
sarcomereLength                        = 2*sarcomereProperties.actinLength;




%One attached crossbridge
kamMin = calcActinMyosinLoadPathStiffness(...
                      1, ...
                      lengthMLineToMyosinMeanAttachmentPoint, ...
                      sarcomereLength,...
                      sarcomereProperties);

%All 20% of the possible attached crossbridges, which is the maximum
%according to Howard 1997
kamMax = calcActinMyosinLoadPathStiffness(...
                      numberOfAttachedCrossbridges, ...
                      lengthMLineToMyosinMeanAttachmentPoint, ...
                      sarcomereLength,...
                      sarcomereProperties);


%Titin stiffness (from Kellermeyer & Granzier 1997 Fig 4 B)
kt2 = [0.00580,0.0288]; %pN/nm at 2um
kt4 = [0.05050,0.0928]; %pN/nm at 4um

kt2T = kt2.*6; %There are 6 titins per half sarcomere
kt4T = kt4.*6; %There are 6 titins per half sarcomere

kt2Ta = kt2T.*ratio_kTitinActive_kTitinPassive;
kt4Ta = kt4T.*ratio_kTitinActive_kTitinPassive;

%%
% Plot the data
%%
numberOfHorizontalPlotColumns = 1;
numberOfVerticalPlotRows      = 1;
flag_fixedPlotWidth           = 1;
plotWidth                     = 4;
plotHeight                    = 7;
pageWidth = plotWidth+5;
pageHeight = plotHeight+5;
plotHorizMarginCm = 2;
plotVertMarginCm = 2;

flag_usingOctave              = 0;
plotConfigGeneric;

fig_loadPathStiffness = figure;

subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

xLineColor = [0,0,0];
tLineColor = [0,0,1];
taLineColor = [1,0,0];


dy = 0.1;

yTop    = max(kamMax) + dy*(20*max(kamMax));
yBottom = 10.^(round(min(log10(kt2T))-1));

yTextA= max(kamMax) + dy*(9*max(kamMax)); 
yTextB= max(kamMax) + dy*(1*max(kamMax));

yTextC = max(kt4T) + dy*(9*max(kt4T));
yTextD = max(kt4T) + dy*(1*max(kt4T));

yTextE = max(kt4Ta) + dy*(9*max(kt4Ta));
yTextF = max(kt4Ta) + dy*(1*max(kt4Ta));


x0 = 1;
x1 = 0.125;
y0 = kamMin(1,1);
y1 = kamMin(1,2);
semilogy([x0,x0],[y0,y1],'Color',xLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',xLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',xLineColor,'LineWidth',1);
hold on;

text(x0,yTextB,'1 XB','HorizontalAlignment','center',...
                  'VerticalAlignment','bottom');
hold on;


text(x0+x1,y0,sprintf('%1.1e',y0),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;
text(x0+x1,y1,sprintf('%1.1e',y1),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;


x0 = 2;
y0 = kamMax(1,1);
y1 = kamMax(1,2);
semilogy([x0,x0],[y0,y1],'Color',xLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',xLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',xLineColor,'LineWidth',1);
hold on;


text(x0,yTextB,sprintf('%1.1f XB',numberOfAttachedCrossbridges),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
hold on;
text(x0+x1,y0,sprintf('%1.1e',y0),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;
text(x0+x1,y1,sprintf('%1.1e',y1),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;


x0 = 3;
y0 = kt2T(1,1);
y1 = kt2T(1,2);
semilogy([x0,x0],[y0,y1],'Color',tLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',tLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',tLineColor,'LineWidth',1);
hold on;

text(x0,yTextD,'2$\mu$m ','HorizontalAlignment','center',...
                          'VerticalAlignment','bottom',...
                          'Color',[0,0,1]);
hold on;
text(x0+x1,y0,sprintf('%1.1e',y0),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;
text(x0+x1,y1,sprintf('%1.1e',y1),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;


x0 = 4;
y0 = kt4T(1,1);
y1 = kt4T(1,2);
semilogy([x0,x0],[y0,y1],'Color',tLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',tLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',tLineColor,'LineWidth',1);
hold on;
text(x0,yTextD,'4$\mu$m ','HorizontalAlignment','center', ...
                         'VerticalAlignment','bottom',...
                          'Color',[0,0,1]);
hold on;
text(x0+x1,y0,sprintf('%1.1e',y0),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;
text(x0+x1,y1,sprintf('%1.1e',y1),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
hold on;



x0 = 5;
y0 = kt2Ta(1,1);
y1 = kt2Ta(1,2);
semilogy([x0,x0],[y0,y1],'Color',taLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',taLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',taLineColor,'LineWidth',1);
hold on;

text(x0,yTextF,'2$\mu$m ','HorizontalAlignment','center',...
                          'VerticalAlignment','bottom',...
                          'Color',taLineColor);
hold on;
text(x0+x1,y0,sprintf('%1.1e',y0),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle',...
     'Color',[0,0,0]);
hold on;
text(x0+x1,y1,sprintf('%1.1e',y1),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle',...
     'Color',[0,0,0]);
hold on;


x0 = 6;
y0 = kt4Ta(1,1);
y1 = kt4Ta(1,2);
semilogy([x0,x0],[y0,y1],'Color',taLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',taLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',taLineColor,'LineWidth',1);
hold on;
text(x0,yTextF,'4$\mu$m ','HorizontalAlignment','center', ...
                         'VerticalAlignment','bottom',...
                          'Color',taLineColor);

hold on;
text(x0+x1,y0,sprintf('%1.1e',y0),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle',...
     'Color',[0,0,0]);
hold on;
text(x0+x1,y1,sprintf('%1.1e',y1),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle',...
     'Color',[0,0,0]);
hold on;

text(1.75,yTextA,'Actin-myosin (AM)','Color',[0,0,0],...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
hold on;

text(3.5,yTextC,'Passive Titin (TP)','Color',[0,0,1],...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
hold on;

text(5.5,yTextE,'Active Titin (TA)','Color',taLineColor,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
hold on;


xlim([0.5,6.5]);
xticks([1,2,3,4,5,6]);
xticklabels({'AM: Low', 'AM: High', 'TP: Low', 'TP: High', 'TA: Low', 'TA: High'});
xtickangle(45);
xlabel('Load Paths');

ylabel('Stiffness (pN/nm)');
yticks([0.01,0.1,1,10,100]);


box off;
ylim([yBottom,yTop]);
text(-0.75,yTop,...
  'Stiffness comparison: actin-myosin \& titin','FontSize',12,...
  'VerticalAlignment','bottom');
%title('Stiffness comparision','HorizontalAlignment','left');
hold on;


print('-dpdf', fullfile(projectFolders.output_plots,...
                'fig_Pub_StiffnessActinMyosinVsTitin.pdf'));
