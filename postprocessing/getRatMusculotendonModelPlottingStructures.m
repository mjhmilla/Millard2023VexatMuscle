function [dataToPlot, dataIndexes, plotSettings] = ...
    getRatMusculotendonModelPlottingStructures(muscleName)

%
% Layout
%
titleMuscleName = '';
switch muscleName
    case 'SOL'
        titleMuscleName = 'Rat Soleus';
    case 'EDL'
        titleMuscleName = 'Rat EDL';        
    otherwise assert(0,'Error muscleName not recognized');
end

plotSettings(2) = ...
    struct('row',0,...
           'col',0,...
           'xlim',[],...
           'ylim',[],...
           'xlabel','',...
           'ylabel','',...
           'title','',...
           'legendLocation','',...
           'xticks',[],...
           'yticks',[]);
idx=1;
plotSettings(idx).row = 1;
plotSettings(idx).col = 1;
plotSettings(idx).xlim = [1.28,4.44];
plotSettings(idx).ylim = [0,1.6];
plotSettings(idx).xlabel = 'Length ($$\mu$$m)';
plotSettings(idx).ylabel = 'Norm. Force ($$f/f_o^M$$)';
plotSettings(idx).title = {[titleMuscleName,' $$f^L(\ell^M)$$']};
plotSettings(idx).legendLocation = 'NorthWest';
plotSettings(idx).xticks = [1.28,1.81,2.53];
plotSettings(idx).yticks = [0.00,0.54,1.00];



idx=idx+1;
plotSettings(idx).row = 1;
plotSettings(idx).col = 2;
plotSettings(idx).xlim = [];
plotSettings(idx).ylim = [0,1.6];
plotSettings(idx).xlabel = 'Velocity ($$\ell^M/\ell_o^M$$)';
plotSettings(idx).ylabel = 'Norm. Force ($$f/f_o^M$$)';
plotSettings(idx).title = {[titleMuscleName,' $$f^{V}(v^M/v^M_{max})$$']};
plotSettings(idx).legendLocation = 'SouthEast';
plotSettings(idx).xticks = [];
plotSettings(idx).yticks = [0,1,1.24,1.44];



switch muscleName
    case 'SOL'
        plotSettings(idx).xlim = [-1.03,1.03];
        plotSettings(idx).xticks = [-1.02,0,1.02];
    case 'EDL'
        plotSettings(idx).xlim = [-2.26,2.26];
        plotSettings(idx).xticks = [-2.25,0,2.25];

    otherwise assert(0,'Error: muscleName not found');
end

%
% Index series
%

dataIndexes.SW1982_fl           = nan;
dataIndexes.SW1982_fpe          = nan;
dataIndexes.ZHGL1995_fl         = nan;
dataIndexes.TRSS2017_fl         = nan;
dataIndexes.TRSS2017_fpe        = nan;
dataIndexes.model_fl            = nan;
dataIndexes.model_titinPassive  = nan;
dataIndexes.model_titinActive   = nan;
%dataIndexes.model_titinActiveKmax   = nan;
dataIndexes.model_fv            = nan;
dataIndexes.model_fpe           = nan;

%
% Data series
%

% 1. Tomalka et al. (fl & fpe)
% 2. Zuurbier et al. (fl)
% 3. Stephenson and Williams (fl & fpe)
dataToPlot(10)=...
    struct('x',[],...
           'y',[],...
           'row',0,...
           'col',0,...
           'type','Exp',...
           'LineColor',[0,0,0],...
           'Mark','-',...
           'MarkerFaceColor',[0,0,0],...
           'MarkerEdgeColor',[0,0,0],...
           'MarkerSize',5,...
           'DisplayName','',...
           'HandleVisibility','off',...
           'enablePlot',0);

%1. SW1992 fl
%2. SW1992 fpe
%3. ZHGL1995 fl
%4. TRSS2017 fl
%5. TRSS2017 fpe
%6. fl (model)
%7. fpe (model)
%8. ft passive (model)
%9. fv (model)

modelColors = getPaulTolColourSchemes('vibrant');
colorTRSS2017   = [1,1,1].*0;
colorZHGL1995   = [1,1,1].*0.75;
colorSW1992     = [1,1,1].*0.5;

%SW1982 fl
idx=1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Exp';
dataToPlot(idx).LineColor=colorSW1992;
dataToPlot(idx).Mark='o';
dataToPlot(idx).MarkerFaceColor=[1,1,1];
dataToPlot(idx).MarkerEdgeColor=colorSW1992;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='SW1982 (F)';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.SW1982_fl           = idx;


%SW1992 fpe
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Exp';
dataToPlot(idx).LineColor=colorSW1992;
dataToPlot(idx).Mark='o';
dataToPlot(idx).MarkerFaceColor=[1,1,1];
dataToPlot(idx).MarkerEdgeColor=colorSW1992;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='SW1982 (F)';
dataToPlot(idx).HandleVisibility = 'off';
dataToPlot(idx).enablePlot = 1;

dataIndexes.SW1982_fpe          = idx;


%ZHGL1995 fl
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Exp';
dataToPlot(idx).LineColor=colorZHGL1995;
dataToPlot(idx).Mark='o';
dataToPlot(idx).MarkerFaceColor=colorZHGL1995;
dataToPlot(idx).MarkerEdgeColor=colorZHGL1995;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='ZHGL1995 (B)';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.ZHGL1995_fl         = idx;


%TRSS2017 fl
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Exp';
dataToPlot(idx).LineColor=colorTRSS2017;
dataToPlot(idx).Mark='.';
dataToPlot(idx).MarkerFaceColor=colorTRSS2017;
dataToPlot(idx).MarkerEdgeColor=colorTRSS2017;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='TRSS2017 (F)';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.TRSS2017_fl         = idx;



%TRSS2017 fpe
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Exp';
dataToPlot(idx).LineColor=colorTRSS2017;
dataToPlot(idx).Mark='.';
dataToPlot(idx).MarkerFaceColor=colorTRSS2017;
dataToPlot(idx).MarkerEdgeColor=colorTRSS2017;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='TRSS2017 (F)';
dataToPlot(idx).HandleVisibility = 'off';
dataToPlot(idx).enablePlot = 1;

dataIndexes.TRSS2017_fpe        = idx;


%fl (model)
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Mdl';
dataToPlot(idx).LineColor=modelColors.blue;
dataToPlot(idx).Mark='-';
dataToPlot(idx).MarkerFaceColor=modelColors.blue;
dataToPlot(idx).MarkerEdgeColor=modelColors.blue;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='$$f^L(\ell^M)$$';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.model_fl            = idx;


%ft passive (model)
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Mdl';
dataToPlot(idx).LineColor=modelColors.cyan;
dataToPlot(idx).Mark='-';
dataToPlot(idx).MarkerFaceColor=modelColors.cyan;
dataToPlot(idx).MarkerEdgeColor=modelColors.cyan;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='$$f^{Ti}(\ell^M)$$ (passive)';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.model_titinPassive  = idx;


%ft active (model: titin-actin bond at N2A)
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Mdl';
dataToPlot(idx).LineColor=modelColors.magenta;
dataToPlot(idx).Mark='-';
dataToPlot(idx).MarkerFaceColor=modelColors.magenta;
dataToPlot(idx).MarkerEdgeColor=modelColors.magenta;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='$$f^{Ti}(\ell^M)$$ (active)';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.model_titinActive   = idx;

% %ft active (model: titin-actin bond at PEVK-IgD)
% idx=idx+1;
% dataToPlot(idx).row=1;
% dataToPlot(idx).col=1;
% dataToPlot(idx).x = [];
% dataToPlot(idx).y = [];
% dataToPlot(idx).type = 'Mdl';
% dataToPlot(idx).LineColor=modelColors.orange;
% dataToPlot(idx).Mark='-';
% dataToPlot(idx).MarkerFaceColor=modelColors.orange;
% dataToPlot(idx).MarkerEdgeColor=modelColors.orange;
% dataToPlot(idx).MarkerSize=5;
% dataToPlot(idx).DisplayName='$$f^{Ti}(\ell^M)$$ (active: $$k_{max}$$)';
% dataToPlot(idx).HandleVisibility = 'on';
% 
% dataIndexes.model_titinActiveKmax   = idx;


%fv (model)
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=2;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Mdl';
dataToPlot(idx).LineColor=modelColors.blue;
dataToPlot(idx).Mark='-';
dataToPlot(idx).MarkerFaceColor=modelColors.blue;
dataToPlot(idx).MarkerEdgeColor=modelColors.blue;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='$$f^V(v^M)$$';
dataToPlot(idx).HandleVisibility = 'on';
dataToPlot(idx).enablePlot = 1;

dataIndexes.model_fv            = idx;

%fpe (model)
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Mdl';
dataToPlot(idx).LineColor=modelColors.cyan;
dataToPlot(idx).Mark='--';
dataToPlot(idx).MarkerFaceColor=modelColors.cyan;
dataToPlot(idx).MarkerEdgeColor=modelColors.cyan;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='$$f^{PE}(\ell^M)$$';
dataToPlot(idx).enablePlot = 0;

dataIndexes.model_fpe           = idx;
