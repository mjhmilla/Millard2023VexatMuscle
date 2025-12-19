clc; 
close all;
clear all;

%%
%
% Before running this script, make sure that
% 1. Open main_CreateRatMuscleModel.m and set
%    flag_OuterLoopMode=1;
%    on line  8
% 2. Open main_TomalkaRodeSchumacherSiebert2017.m and set
%    flag_OuterLoopMode=1;
%    on line  7
%
% Once the script is run, you will find the following outputs:
%
% 1. Figure 2 from the paper
%       output/plots/MuscleCurves/fig_Pub_RatMuscleCurves_TRSS2017_0.pdf
%
%
% 2. Figure 4 and 5
%
%       Figure 4 is in the 2nd row while Figure 5A,B,C are in the 3rd row:
%       output/structs/plots/TomalkaRodeSchumacherSiebert2017/
%           fig_Sim_TRSS2017_123_Fl_Fv_QToF.pdf
%
%       Figure 4 is in the 2nd row (yes, a duplicate) while Figure 5D,D,F 
%       are in the 3rd row:
%       output/structs/plots/TomalkaRodeSchumacherSiebert2017/
%           fig_Sim_TRSS2017_123_Fl_Fv_QToF_i.pdf
%
% 3. Rat EDL VEXAT model parameters (before fitting)
%       output/structs/FittedModels/ratTRSS2017EDLFibrilActiveTitin_0.mat
%       output/structs/FittedModels/ratTRSS2017EDLFibrilActiveTitin_1.mat
%       output/structs/FittedModels/ratTRSS2017EDLFibrilActiveTitin_2.mat
%    *Note: These three structs are identical.
%   
% 4. Rat EDL VEXAT model parameters (after fitting)
%
%   One single fitted Q value for all 3 trials Q
%       output/structs/FittedModels
%           /ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF.mat
%
%   Individually fit Q values (Q_1, Q_2, Q_3)
%       output/structs/FittedModels
%           /ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF_i.mat
% 
% 5. Simulation output
%
%   A struct with time series recordings from many of the VEXAT model's
%   internal variables using a single value of Q that's best fit to 
%   the 3 mean fiber responses in Tomalka et al. 2017
%       output/structs/TomalkaRodeSchumacherSiebert2017
%           /benchRecordVexat_TRSS2017_fitted_123_Fl_Fv_QToF.mat
%
%   A struct with "..." with values of Q that are individually fit 
%   (Q_1, Q_2, Q_3) to the three mean fiber responses trials of 
%   Tomalka et al. 2017.
%       output/structs/TomalkaRodeSchumacherSiebert2017
%           /benchRecordVexat_TRSS2017_fitted_123_Fl_Fv_QToF_i.mat
%
% 6. Fitting information
%   
%   The information that appears in Table 1 when a single Q value is fit
%   to all three trials can be found here:
%       fittingLog__123_Fl_Fv_QToF.txt
%       fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF.mat
%
%   The information that appears in Table 1 when each of the three trials 
%   have individually fit Q values can be found here:
%       fittingLog__123_Fl_Fv_QToF_i.txt   
%       fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF_i.mat
%%

main_CreateRatMuscleModel;

simConfigInput.runFitting              = 1; 
simConfigInput.generatePlots           = 1;
simConfigInput.fitToIndividualTrials   = 1; %Individually fitted Q's
simConfigInput.manuallySetTimeConstant = 0;

main_TomalkaRodeSchumacherSiebert2017;

clc;
close all;
clear all;

simConfigInput.runFitting              = 1; 
simConfigInput.generatePlots           = 1;
simConfigInput.fitToIndividualTrials   = 0;  %One Q for all trials
simConfigInput.manuallySetTimeConstant = 0;

main_TomalkaRodeSchumacherSiebert2017;

