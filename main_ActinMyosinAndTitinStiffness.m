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


% Huxley & Simmons 1981: number crossbridges per half-sarcomere
numCrossbridgesHalfSarcomere = 2*(700/14.3); %pair of xbridges every 14.3 nm for 700 nm

%One attached crossbridge
kam1 = calcActinMyosinLoadPathStiffness(...
                      1, ...
                      lengthMLineToMyosinMeanAttachmentPoint, ...
                      sarcomereLength,...
                      sarcomereProperties);

%All 20% of the possible attached crossbridges, which is the maximum
%according to Howard 1997
kam20p = calcActinMyosinLoadPathStiffness(...
                      0.20*numCrossbridgesHalfSarcomere, ...
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

yTop = max(kam20p)*1.1;
yBottom=min(kt2T);


yTextA=55;%0.03;
yTextB=37.5;%0.02;
yTextC = 1.75;
yTextD = 1.1;
yTextE = 3.5;
yTextF = 2.1;


x0 = 1;
x1 = 0.125;
y0 = kam1(1,1);
y1 = kam1(1,2);
semilogy([x0,x0],[y0,y1],'Color',xLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',xLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',xLineColor,'LineWidth',1);
hold on;

text(x0,yTextB,'1 XB','HorizontalAlignment','center',...
                  'VerticalAlignment','top');
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
y0 = kam20p(1,1);
y1 = kam20p(1,2);
semilogy([x0,x0],[y0,y1],'Color',xLineColor,'LineWidth',2);
hold on;
semilogy([x0-x1,x0+x1],[y0,y0],'Color',xLineColor,'LineWidth',1);
hold on;
semilogy([x0-x1,x0+x1],[y1,y1],'Color',xLineColor,'LineWidth',1);
hold on;


text(x0,yTextB,'19.6 XB','HorizontalAlignment','center',...
                        'VerticalAlignment','top');
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
                          'VerticalAlignment','top',...
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
                         'VerticalAlignment','top',...
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
                          'VerticalAlignment','top',...
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
                         'VerticalAlignment','top',...
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

text(1.5,yTextA,'Actin-myosin (AM)','Color',[0,0,0],...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
hold on;

text(3.5,yTextC,'Passive Titin (TP)','Color',[0,0,1],...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
hold on;

text(5.5,yTextE,'Active Titin (TA)','Color',taLineColor,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');
hold on;


xlim([0.5,6.5]);
xticks([1,2,3,4,5,6]);
xticklabels({'AM: Low', 'AM: High', 'TP: Low', 'TP: High', 'TA: Low', 'TA: High'});
xtickangle(45);
xlabel('Load Paths');

ylabel('Stiffness (pN/nm)');
yticks([0.01,0.1,1,10]);


box off;

text(-0.75,75,'Stiffness comparison: actin-myosin \& titin','FontSize',12);
%title('Stiffness comparision','HorizontalAlignment','left');
hold on;
ylim([0.011,55]);

print('-dpdf', fullfile(projectFolders.output_plots,...
                'fig_Pub_StiffnessActinMyosinVsTitin.pdf'));
