function [R_Index Art_Index] = DetectR(Data,SamplingRate,HR_Range)
% This function detects the R-waves in ECG signal
%
% Inputs
% Data             : ECG time series data (a vector)
% SamplingRate     : Sampling rate (Number of samples per second)
% HR_Range         : Normal range of heart rate : [60 100] adult deafult
%
% Output
% R_Index          : Index of each detected R wave
% Art_Index        : Index of R_Wave's that are artifacts 
%
% Copyright (c) 2018, Farhad Kaffashi & Kenneth A. Loparo
% All rights reserved.
% Distribution without permission is prohibited.
%
% Please contact Farhad Kaffashi <farhad@case.edu> or Kenneth A. Loparo <kenneth.loparo@case.edu>
% if you have any questions.
