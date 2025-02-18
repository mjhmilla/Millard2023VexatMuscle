function [dataToPlot, dataIndexes, plotSettings] = getRatSoleusModelPlottingStructures()

%
% Layout
%

plotSettings(2) = ...
    struct('row',0,...
           'col',0,...
           'xlabel','',...
           'ylabel','',...
           'xlim',[],...
           'ylim',[],...
           'title','');
idx=1;
plotSettings(idx).row = 1;
plotSettings(idx).col = 1;
plotSettings(idx).xlim = [1.3,4.5];
plotSettings(idx).ylim = [0,1.6];
plotSettings(idx).xlabel = 'Length ($$\mu$$m)';
plotSettings(idx).ylabel = 'Norm. Force ($$f/f_o^M$$)';
plotSettings(idx).title = {'Rat Soleus $$f^L(\ell^M)$$'};
plotSettings(idx).legendLocation = 'NorthWest';

% idx=idx+1;
% plotSettings(idx).row = 1;
% plotSettings(idx).col = 2;
% plotSettings(idx).xlim = [1.3,4.5];
% plotSettings(idx).ylim = [0,1.1];
% plotSettings(idx).xlabel = 'Length ($$\mu$$m)';
% plotSettings(idx).ylabel = 'Norm. Force ($$f/f_o^M$$)';
% plotSettings(idx).title = {'Rat Soleus $$f^{PE}(\ell^M)$$'};

idx=idx+1;
plotSettings(idx).row = 1;
plotSettings(idx).col = 2;
plotSettings(idx).xlim = [-1,1];
plotSettings(idx).ylim = [0,1.6];
plotSettings(idx).xlabel = 'Velocity ($$\ell^M/\ell_o^M$$)';
plotSettings(idx).ylabel = 'Norm. Force ($$f/f_o^M$$)';
plotSettings(idx).title = {'Rat Soleus $$f^{V}(v^M/v^M_{max})$$'};
plotSettings(idx).legendLocation = 'SouthEast';

%
% Index series
%

dataIndexes.SW1982_fl           = 1;
dataIndexes.SW1982_fpe          = 2;
dataIndexes.ZHGL1995_fl         = 3;
dataIndexes.TRSS2017_fl         = 4;
dataIndexes.TRSS2017_fpe        = 5;
dataIndexes.model_fl            = 6;
dataIndexes.model_titinPassive  = 7;
dataIndexes.model_titinActive   = 8;
dataIndexes.model_fv            = 9;
dataIndexes.model_fpe           = 10;

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
           'HandleVisibility','off');

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

%SW1992 fl
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
dataToPlot(idx).DisplayName='SW1992 (F)';
dataToPlot(idx).HandleVisibility = 'on';

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
dataToPlot(idx).DisplayName='SW1992 (F)';
dataToPlot(idx).HandleVisibility = 'off';

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


%ft active (model)
idx=idx+1;
dataToPlot(idx).row=1;
dataToPlot(idx).col=1;
dataToPlot(idx).x = [];
dataToPlot(idx).y = [];
dataToPlot(idx).type = 'Mdl';
dataToPlot(idx).LineColor=modelColors.red;
dataToPlot(idx).Mark='-';
dataToPlot(idx).MarkerFaceColor=modelColors.red;
dataToPlot(idx).MarkerEdgeColor=modelColors.red;
dataToPlot(idx).MarkerSize=5;
dataToPlot(idx).DisplayName='$$f^{Ti}(\ell^M)$$ (active)';
dataToPlot(idx).HandleVisibility = 'on';

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
dataToPlot(idx).DisplayName='$$f^V(v^M/v_o^M)$$';
dataToPlot(idx).HandleVisibility = 'on';


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
