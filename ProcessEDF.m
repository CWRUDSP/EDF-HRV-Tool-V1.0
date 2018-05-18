%--------------------------------------------------------------------------
% @license
% Copyright 2018 IDAC Signals Team, Case Western Reserve University 
%
% Lincensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public 
% you may not use this file except in compliance with the License.
%
% Unless otherwise separately undertaken by the Licensor, to the extent possible, 
% the Licensor offers the Licensed Material as-is and as-available, and makes no representations 
% or warranties of any kind concerning the Licensed Material, whether express, implied, statutory, or other. 
% This includes, without limitation, warranties of title, merchantability, fitness for a particular purpose, 
% non-infringement, absence of latent or other defects, accuracy, or the presence or absence of errors, 
% whether or not known or discoverable. 
% Where disclaimers of warranties are not allowed in full or in part, this disclaimer may not apply to You.
%
% To the extent possible, in no event will the Licensor be liable to You on any legal theory 
% (including, without limitation, negligence) or otherwise for any direct, special, indirect, incidental, 
% consequential, punitive, exemplary, or other losses, costs, expenses, or damages arising out of 
% this Public License or use of the Licensed Material, even if the Licensor has been advised of 
% the possibility of such losses, costs, expenses, or damages. 
% Where a limitation of liability is not allowed in full or in part, this limitation may not apply to You.
%
% The disclaimer of warranties and limitation of liability provided above shall be interpreted in a manner that, 
% to the extent possible, most closely approximates an absolute disclaimer and waiver of all liability.
%
% Developed by the IDAC Signals Team at Case Western Reserve University 
% with support from the National Institute of Neurological Disorders and Stroke (NINDS) 
%     under Grant NIH/NINDS U01-NS090405 and NIH/NINDS U01-NS090408.
%              Farhad Kaffashi
%--------------------------------------------------------------------------
function [R_Index Art_Index SelectedCh SamplingRate]=ProcessEDF(FileName,FileInfo,ChInfo,HR_Range)
%------------------------------------------------------------------
%% find ECG or EKG Ch
TempCh = [];
Index  = [];
for i=1:FileInfo.SignalNumbers
    TempCh{i} = ChInfo.Labels(i,:);
    Temp1 = find(TempCh{i}>=65 & TempCh{i}<=90);
    TempCh{i}(Temp1) = TempCh{i}(Temp1) + 32;
    Temp1 = strfind(TempCh{i},'ecg');
    if ~isempty(Temp1)
        Index = [Index i];
    end
    Temp1 = strfind(TempCh{i},'ekg');
    if ~isempty(Temp1)
        Index = [Index i];
    end
    
    Temp1 = strfind(TempCh{i},'x1');
    if ~isempty(Temp1)
        Index = [Index i];
    end
    Temp1 = strfind(TempCh{i},'x2');
    if ~isempty(Temp1)
        Index = [Index i];
    end
end

if ~isempty(Index)
    
    SamplingRate = fix(ChInfo.nr(Index(1))/FileInfo.DataRecordDuration);
    
    if SamplingRate < 500
        %------------------------------------------------------------------
        %% Read the data file
        
        SelectedCh = Index;
        Index = [0 cumsum(ChInfo.nr)'];
        % Read the data for processing
        
        tic
        Fid = fopen(FileName,'r');
        
        fseek(Fid,FileInfo.HeaderNumBytes ,'bof');
        TotalEKG=zeros(length(SelectedCh),FileInfo.DataRecordDuration*FileInfo.NumberDataRecord*SamplingRate);

        % Read one minute data at a time
        NumEpoch = fix(60/FileInfo.DataRecordDuration);
        %%
        
        Temp = fread(Fid,[sum(ChInfo.nr) NumEpoch],'int16');
        ShiftSample = 0;
        
        while ~isempty(Temp)
            
            Counter = 0;
            for i=SelectedCh
                Counter = Counter + 1;
                Temp1 = Temp([1:ChInfo.nr(i)]+Index(i),:);
                Temp1 = Temp1(:)';
                TotalEKG(Counter,[1:length(Temp1)]+ShiftSample) = Temp1;
            end
            
            ShiftSample = ShiftSample + length(Temp1);
            Temp = fread(Fid,[sum(ChInfo.nr) NumEpoch],'int16');
        end
        
        
        
        
        fclose(Fid);
        
        
        
        
        %% Detect the R Waves and Artifact
        
        R_Index   = [];
        Art_Index = [];
        for i = 1:length(SelectedCh)
            
            [Temp1 Temp2] = DetectR_Art(TotalEKG(i,:),SamplingRate,HR_Range);
            R_Index{i}   = Temp1;
            Art_Index{i} = Temp2;
            
        end
        toc
        
    end
end
