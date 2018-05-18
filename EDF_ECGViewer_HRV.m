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
function varargout = EDF_ECGViewer_HRV(varargin)
% EDF_ECGVIEWER_HRV MATLAB code for EDF_ECGViewer_HRV.fig
%      EDF_ECGVIEWER_HRV, by itself, creates a new EDF_ECGVIEWER_HRV or raises the existing
%      singleton*.
%
%      H = EDF_ECGVIEWER_HRV returns the handle to a new EDF_ECGVIEWER_HRV or the handle to
%      the existing singleton*.
%
%      EDF_ECGVIEWER_HRV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDF_ECGVIEWER_HRV.M with the given input arguments.
%
%      SAFEBEATDATAVI EWER_ANN_MAC('Property','Value',...) creates a new EDF_ECGVIEWER_HRV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EDF_ECGViewer_HRV_OpeningFcn gets calclcled.  An
%      unrecognized property name or detectinvalid value makes property application
%      stop.  All inputs are passed to EDF_ECGViewer_HRV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help EDF_ECGViewer_HRV

% Last Modified by GUIDE v2.5 30-Mar-2018 11:51:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EDF_ECGViewer_HRV_OpeningFcn, ...
    'gui_OutputFcn',  @EDF_ECGViewer_HRV_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EDF_ECGViewer_HRV is made visible.
function EDF_ECGViewer_HRV_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EDF_ECGViewer_HRV (see VARARGIN)

% Choose default command line output for EDF_ECGViewer_HRV
handles.output = hObject;

handles.FlagAnn = 0;
handles.FlagCliper = 0;
handles.FlagSaveFile = 0;

handles.OldPlotMode.Trace   = 4; % Ch1
handles.OldPlotMode.Signals = 1; % Ch1 & Ch2


% Sensitivity Sel for each Mode
% 1 for regular mode
% 2 for trace mode
% 3 default value 3 for both mode
handles.SenSel = [8 8];

handles.SelectedFolder = 1;

handles.UserAnn = [];

handles.AxesPosition = get(handles.axes1,'position');
handles.PlotMode = 1 ;

handles.MinMaxValue = [500 1500];

handles.AnnType = {'Artifact',...
    'Signal Break',...
    'Lost ECG Lead Start',...
    'Lost ECG Lead End' ...
    'PAC',...
    'Aberrant PAC',...
    'Atrial Fibrillation',...
    'Atrial Flutter',...
    'Start Atrial Fibrillation',...
    'End Atrial Fibrillation', ...
    'SVT',...
    'PVC',...
    'Ventricular Bigeminy',...
    'Ventricular Trigeminy',...
    'Ventricular Triplet',...
    'Ventricular Quadrigeminy',...
    'Ventricular Couplet', ...
    'NSVT',...
    'Ventricular Tachycardia',...
    'IVCD',...
    'First Degree A-V Block',...
    'Second Degree A-V Block',...
    'Third Degree A-V Block',...
    'Sinus Pause',...
    'Sinus Arrest/ Ventricular Asystole',...
    'Sinus Tachycardia',...
    'Sinus Bradycardia',...
    'Idioventricular Rhythm',...
    'Motion Artifact Start',...
    'Motion Artifact End'};

handles.Ann = [];

set(handles.ListBoxAnnType,'string',handles.AnnType);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EDF_ECGViewer_HRV wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.axes1,'xTickLabel','','yTickLabel','')
set(handles.axes2,'xTickLabel','','yTickLabel','')
set(handles.axes3,'xTickLabel','','yTickLabel','')
set(handles.axes4,'xTickLabel','','yTickLabel','')
set(handles.axes5,'xTickLabel','','yTickLabel','')



% --- Outputs from this function are returned to the command line.
function varargout = EDF_ECGViewer_HRV_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.FlagSaveFile == 1
    warndlg('Save The Annotation Before Changing the File')
else
    handles.Path = uigetdir;
    
    if ispc
        handles.Path = [handles.Path '\'];
    else
        handles.Path = [handles.Path '/'];
    end
    
    % find all folders under the given Main folder
    Temp = dir([handles.Path]);
    
    if ~(length(Temp)==0)
        % detect all the paths
        handles.FolderNames = [];
        TempNames = [];
        Counter = 0;
        for i=1:length(Temp)
            if Temp(i).isdir
                if length(Temp(i).name)>2
                    Counter = Counter + 1;
                    handles.FolderNames{Counter} = Temp(i).name;
                    TempNames{Counter} = [num2str(Counter) '. ' Temp(i).name];
                end
            end
        end
    end
    set(handles.ListBoxFolderNames,'string',TempNames,'value',1);
    guidata(hObject,handles)
end


% --- Executes on selection change in ListBoxFolderNames.
function handles = ListBoxFolderNames_Callback(hObject, eventdata, handles)
% hObject    handle to ListBoxFolderNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for CounterFile = get(hObject,'value')
    %for CounterFile =  266 : length(get(hObject,'string'))
    
    set(hObject,'value',CounterFile);
    
    
    %%
    
%     
%     if (handles.FlagSaveFile == 1) | ~(get(handles.CheckBoxArtifact,'value'))
%         %warndlg('Save The Annotation Before Changing the File')
%         MinMaxValue=handles.MinMaxValue;
%         UserAnn = handles.UserAnn;
%         FlagArtifact = get(handles.CheckBoxArtifact,'value');
%         
%         
%         if ispc
%             Path = [handles.Path handles.FolderNames{handles.SelectedFolder} '\'];
%         else
%             Path = [handles.Path handles.FolderNames{handles.SelectedFolder} '/'];
%         end
%         
%         FileName  = [Path handles.FileNames{get(handles.ListBoxFileNames,'value')}([1:end-4]) '_Ann'];
%         
%         save(FileName,'MinMaxValue','UserAnn','FlagArtifact');
%         
%         handles.FlagSaveFile = 0;
%         clear FlagArtifact;
%         handles.FlagSaveFile
%         set(handles.CheckBoxArtifact,'value',1);
%     end
%     
    
    %%
    
    handles.UserAnn = [];
    Temp = get(hObject,'value');
    handles.SelectedFolder = Temp;
    
    
    % Get the file names
    
    if ispc
        Path = [handles.Path handles.FolderNames{handles.SelectedFolder} '\'];
    else
        Path = [handles.Path handles.FolderNames{handles.SelectedFolder} '/'];
    end
    
    % find all EDF Files under the given folder
    Temp = dir([Path '*.EDF']);
    
    if ~(length(Temp)==0)
        % detect all the paths
        handles.FileNames = [];

        TempNames = [];
        for i=1:length(Temp)
            handles.FileNames{i} = Temp(i).name;
            % get the data & time
            FileName=[Path handles.FileNames{i}];
            
            Fid = fopen(FileName);
            
            fseek(Fid,168,'bof');
            Temp1 = char(fread(Fid,[1 16],'uint8'));
            Date=[Temp1(4:5) '/' Temp1(1:2) '/' Temp1(7:8)];
            Time=[Temp1(9:10) ':' Temp1(12:13) ':' Temp1(15:16)];
            
            fclose(Fid);
            
            
            TempNames{i} = [num2str(i) '. ' Date ' ' Time ' ' Temp(i).name];
        end
    end
    set(handles.ListBoxFileNames,'string',TempNames,'value',1);
    guidata(hObject,handles)
    
    
end




function handles=UpdatePlot(handles)

% plot Mode
% 1 : Both Ch1 and Ch2
% 2 : Just Ch1
% 3 : Just Ch2
% 4 : Trace plot ch1
% 5 : Trace plot ch2

Sen = get(handles.PopMenuSen,'string');
Sen = Sen{get(handles.PopMenuSen,'value')};
Sen([-2:0]+end)=[];
Sen = str2num(Sen)/2;

%---------------------
Temp = get(handles.PopMenuWindowTime,'string');
Temp = Temp{get(handles.PopMenuWindowTime,'value')};
Temp([-3:0]+end)=[];
WindowTime = str2num(Temp);

if isempty(WindowTime)
    WindowTime = str2num(get(handles.EditCustomSize,'string'));
end

%----------------------
% Data Load
if get(handles.CheckBoxTrace,'value')
    NumTrace   = str2num(get(handles.EditNumTrace,'string'));
else
    NumTrace   = 1;
end

NumSamples = fix(WindowTime*handles.SamplingRate);
Shift      = round(get(handles.SliderTime,'value'));
Time       = ([1:NumSamples])/handles.SamplingRate;


SkipSample = fix(WindowTime*NumTrace*handles.SamplingRate);
TempCh       = get(handles.ListBoxChSelect,'value');
TempSelected = get(handles.PopMenuHRV_Ch,'value');

Data = DataLoad(handles.Fid,round(Shift/handles.SamplingRate),SkipSample,handles.PlotMode,handles.ECG_Ch(TempCh),handles.FileInfo,handles.ChInfo);

if get(handles.CheckBoxFilter,'value')
    for i=1:size(Data,1)
        Data(i,:)=filtfilt(handles.FilterB,handles.FilterA,Data(i,:));
    end
end
if get(handles.CheckBoxLowPass,'value')
    for i=1:size(Data,1)
        Data(i,:)=filtfilt(handles.FilterB_LowPass,handles.FilterA_LowPass,Data(i,:));
    end
end

if get(handles.CheckBoxNotchFilter,'value')
    for i=1:size(Data,1)
        Data(i,:)=filtfilt(handles.NotchB,handles.NotchA,Data(i,:));
    end
end

%----------------------
%  detect R

Temp = handles.R_Index{TempSelected};

Index=find(Temp>Shift & Temp<(Shift+str2num(get(handles.EditHRTime,'string'))*handles.SamplingRate));

R_Index = Temp(Index)-Shift;

if ~isempty(R_Index)
    Art_Index = handles.Art_Index{TempSelected}-Index(1)+1;
    Art_Index(Art_Index>length(R_Index))=[];
    Art_Index(Art_Index<1)=[];
else
    Art_Index =[];
end


if get(handles.CheckBoxClipping,'value')
    % clipping the saturated signal
    for i=1:size(Data,1)
        Data(i,:) = Data(i,:) - mean(Data(i,:));
    end
    
    Data(Data>Sen)=Sen;
    Data(Data<-Sen)=-Sen;
end

%------------------------

% plot Heart Rate only for plot mode 1 2 3

if 1%handles.PlotMode<4
    axes(handles.axes3);
    if ~isnan(R_Index)
        R_Time = R_Index(1:end-1)/handles.SamplingRate;
        
        %HR = (1./(diff(R_Index)/handles.SamplingRate))*60;
        HR = diff(R_Index)/handles.SamplingRate*1000;
        
        HR(HR>3000)=nan;
        
        hold off
        if ~isempty(Art_Index)
            if Art_Index(1)<=1
                Art_Index(Art_Index<2)=[];
            end
            if ~isempty(Art_Index)
                if Art_Index(end)==length(R_Index)
                    Art_Index(end)=[];
                end
            end
            if ~isempty(Art_Index)
                plot(R_Time(Art_Index),HR(Art_Index),'r*')
                hold on
                plot(R_Time(Art_Index-1),HR(Art_Index-1),'r*')
                
                HR(Art_Index)=nan;
                HR(Art_Index-1)=nan;
            end
        end
        
        
        
        plot(R_Time,HR,'linewidth',1.5,'color','b');%,'marker','*')
        hold on
        
        set(handles.axes3,'fontweight','bold','xticklabel','')
        xlim([0 str2num(get(handles.EditHRTime,'string'))])
        
        if ~isempty(handles.Ann.Time)
            IndexAnn=find(handles.Ann.Time(:,1)>Shift/handles.SamplingRate & handles.Ann.Time(:,1)<(R_Time(end)+Shift/handles.SamplingRate));
            if ~isempty(IndexAnn)
                TimeAnn = handles.Ann.Time(IndexAnn,1)-Shift/handles.SamplingRate;
                hold on
                Temp=get(handles.axes3,'ylim');
                plot([TimeAnn TimeAnn]',Temp,'r')
                hold off
            end
        end
        
        axes(handles.axes5)
        plot(HR(1:end-1),HR(2:end),'*');
        set(handles.axes5,'fontweight','bold')
    else
        cla
    end
    
    
    
end


%%


axes(handles.axes1)
hold off
if WindowTime<30
    TickDiff = 0.2;
elseif WindowTime>70
    TickDiff = 10;
else
    TickDiff = 0.4;
end


handles.PlotAnn(1) = fill([0 0 0 0], [0 0 0 0] ...
    ,[255 100 50]/255,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
handles.FlagAnn = 0;
handles.TextAnn = text(0,0,'','FontWeight','bold','FontSize',12);
hold on

for i=0:NumTrace-1
    plot(Time,Data(1,[1:NumSamples]+i*NumSamples)-mean(Data(1,:))-i*Sen*2,'linewidth',1.5);
    Temp1 = datestr((Time(1)+WindowTime*i)/86400+Shift/handles.SamplingRate/86400+handles.StartTime,'HH:MM:SS');
    text(Time(end)/100,Sen-Sen/8-i*Sen*2,Temp1,'FontWeight','bold','FontSize',12,'fontname','Times New Roman');
end
ylim([-Sen*NumTrace*2+Sen Sen])


Temp = round([Time(1):TickDiff:Time(end)]*5)/5;
set(handles.axes1,'xtick',Temp,'xticklabel','','fontweight','bold',...
    'GridLineStyle','-','yticklabel','','xlim',[0 WindowTime],'xcolor',[225 200 200]/255,'ycolor',[225 200 200]/255 ...
    );

if get(handles.RadioButtonGrid,'value')
    set(handles.axes1,'xgrid','on','ygrid','on');
end


% ---------------------------------------
% plot Ann only for plot mode 1,2,3

if ~isempty(handles.Ann.Time)
    Start = Shift/handles.SamplingRate;
    End   = Start + WindowTime*NumTrace;
    Index  =        find(handles.Ann.Time(:,1)>Start & handles.Ann.Time(:,1)<End)';
    
    if handles.PlotMode<4
        
        
        Index  = [Index find(handles.Ann.Time(:,2)>Start & handles.Ann.Time(:,2)<End)'];
        Index  = [Index find(handles.Ann.Time(:,1)<Start & handles.Ann.Time(:,2)>End)'];
        
        if ~isempty(Index)
            Index = sort(Index);
            Index(diff(Index)==0)=[];
            for i=Index
                AnnSt  = handles.Ann.Time(i,1);
                if AnnSt<Start
                    AnnSt = Start;
                end
                AnnEnd = handles.Ann.Time(i,2);
                if AnnEnd>End
                    AnnEnd = End;
                end
                fill([AnnSt AnnSt AnnEnd AnnEnd]-Start, [-Sen+Sen/100 Sen-Sen/100 Sen-Sen/100 -Sen+Sen/100] ...
                    ,[255 100 50]/255,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
                text(AnnSt-Start,Sen-Sen/4,handles.Ann.Text{i},'FontWeight','bold','FontSize',14);
            end
        end
        
        % find the closest ann
        [Temp1 Index]=min(abs(handles.Ann.Time(:,1) - (Shift+NumSamples/2)/handles.SamplingRate));
        set(handles.ListBoxAnn,'value',Index);
        
        
        
        
        hold off
    else
        % add annotation for trace mode
        Index  =  [Index      find(handles.Ann.Time(:,2)>Start & handles.Ann.Time(:,2)<End)'];
        Index = sort(Index);
        Index(diff(Index)==0)=[];
        
        
        if ~isempty(Index)
            Color = [255 100 50]/255;
            for i=Index
                AnnSt  = handles.Ann.Time(i,1)-Start;
                AnnEnd = handles.Ann.Time(i,2)-Start;
                
                if AnnSt<0
                    AnnSt=0;
                end
                
                Line = fix(AnnSt/WindowTime);
                AnnSt = AnnSt - Line*WindowTime;
                
                AnnEnd = AnnEnd - Line*WindowTime;
                
                if ~(AnnSt ==0)
                    text(AnnSt,Sen-Sen/5-Sen*2*Line,handles.Ann.Text{i},'FontWeight','bold','FontSize',10);
                end
                
                if Color(3)==50/255
                    Color = [130 180 20]/255;
                else
                    Color = [255 100 50]/255;
                end
                
                while AnnEnd>0
                    
                    
                    if AnnEnd<WindowTime
                        fill([AnnSt AnnSt AnnEnd AnnEnd], [-Sen+Sen/100 Sen-Sen/100 Sen-Sen/100 -Sen+Sen/100]-Sen*2*Line ...
                            ,Color,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
                    else
                        fill([AnnSt AnnSt WindowTime WindowTime], [-Sen+Sen/100 Sen-Sen/100 Sen-Sen/100 -Sen+Sen/100]-Sen*2*Line ...
                            ,Color,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
                    end
                    Line = Line + 1;
                    AnnSt = 0;
                    AnnEnd = AnnEnd - WindowTime;
                    
                end
                
                
            end
        end
        
    end
end
%--------------------------------

if handles.PlotMode==1
    
    
    axes(handles.axes2)
    hold off
    handles.PlotAnn(2) = fill([0 0 0 0], [0 0 0 0] ...
        ,[255 100 50]/255,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
    hold on
    ylim([-Sen Sen])
    plot(Time,Data(2,:)-mean(Data(2,:)),'linewidth',1.5);
    ylim([-Sen Sen])
    
    set(handles.axes2,'xtick',Temp,'xticklabel','','fontweight','bold',...
        'GridLineStyle','-','yticklabel','','xlim',[0 WindowTime],'xcolor',[225 200 200]/255,'ycolor',[225 200 200]/255 ...
        );
    
    if get(handles.RadioButtonGrid,'value')
        set(handles.axes2,'xgrid','on','ygrid','on');
    end
    
    if ~isempty(handles.Ann.Time)
        Start = Shift/handles.SamplingRate;
        End   = Start + WindowTime;
        
        Index  =        find(handles.Ann.Time(:,1)>Start & handles.Ann.Time(:,1)<End)' ;
        Index  = [Index find(handles.Ann.Time(:,2)>Start & handles.Ann.Time(:,2)<End)'];
        Index  = [Index find(handles.Ann.Time(:,1)<Start & handles.Ann.Time(:,2)>End)'];
        
        if ~isempty(Index)
            Index = sort(Index);
            Index(diff(Index)==0)=[];
            for i=Index
                AnnSt  = handles.Ann.Time(i,1);
                if AnnSt<Start
                    AnnSt = Start;
                end
                AnnEnd = handles.Ann.Time(i,2);
                if AnnEnd>End
                    AnnEnd = End;
                end
                fill([AnnSt AnnSt AnnEnd AnnEnd]-Start, [-Sen+Sen/100 Sen-Sen/100 Sen-Sen/100 -Sen+Sen/100] ...
                    ,[255 100 50]/255,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
            end
        end
        
        
    end
    
end

if handles.PlotMode<4
    axes(handles.axes1)
    
    R_Index(R_Index<1 | R_Index>NumSamples)=[];
    if ~isempty(Art_Index)
        Art_Index(Art_Index>length(R_Index))=[];
    end
    hold on
    if ~isnan(R_Index)
        plot(Time(R_Index),Sen-Sen/25,'color','r','marker','v','MarkerFaceColor','r')
        hold on
        if ~isempty(Art_Index)
            plot(Time(R_Index(Art_Index)),Sen-Sen/25,'color','k','marker','v','MarkerFaceColor','k')
        end
        HR = fix(60/(mean(diff(R_Index))/handles.SamplingRate));
        text(Time(end)/100,-Sen+Sen/5,[num2str(HR) ' BPM'],'FontWeight','bold','FontSize',12,'fontname','Times New Roman');
    else
        text(Time(end)/100,-Sen+Sen/5,['Lost ECG'],'FontWeight','bold','FontSize',12,'fontname','Times New Roman');
    end
    hold off
    
end

axes(handles.axes3)




% --- Executes on slider movement.
function SliderTime_Callback(hObject, eventdata, handles)
% hObject    handle to SliderTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


handles=UpdatePlot(handles);
guidata(hObject,handles);


% --- Executes on selection change in PopMenuWindowTime.
function PopMenuWindowTime_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenuWindowTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenuWindowTime contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenuWindowTime

Temp = get(handles.PopMenuWindowTime,'string');
Temp = Temp{get(handles.PopMenuWindowTime,'value')};
Temp([-3:0]+end)=[];
WindowTime = str2num(Temp);

if isempty(WindowTime)
    WindowTime = str2num(get(handles.EditCustomSize,'string'));
    set(handles.EditCustomSize,'visible','on');
else
    set(handles.EditCustomSize,'visible','off');
end

max = fix(handles.FileInfo.DataRecordDuration*handles.FileInfo.NumberDataRecord-WindowTime)*handles.SamplingRate;

Temp = get(handles.SliderTime,'value');

if Temp>max
    set(handles.SliderTime,'value',max);
end

set(handles.SliderTime,'max',max);

Temp = [0.2 1]*WindowTime*handles.SamplingRate/max;
set(handles.SliderTime,'SliderStep',Temp);

handles=UpdatePlot(handles);

guidata(hObject,handles);

function EditCustomSize_Callback(hObject, eventdata, handles)
% hObject    handle to EditCustomSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditCustomSize as text
%        str2double(get(hObject,'String')) returns contents of EditCustomSize as a double

PopMenuWindowTime_Callback(handles.PopMenuWindowTime, eventdata, handles)




function Data = DataLoad(Fid,Shift,NumSamples,PlotMode,Ch,FileInfo,ChInfo)

Data = zeros(2,NumSamples);

SamplingRate = fix(ChInfo.nr(Ch(1))/FileInfo.DataRecordDuration);

SkipNumSegment = Shift/FileInfo.DataRecordDuration;

fseek(Fid,SkipNumSegment*2*sum(ChInfo.nr) +FileInfo.HeaderNumBytes,'bof');

NumEpoch = fix(NumSamples/SamplingRate/FileInfo.DataRecordDuration);
%%

Temp = fread(Fid,[sum(ChInfo.nr) NumEpoch],'int16');

Index = [0 cumsum(ChInfo.nr)'];

if ~isempty(Temp)
    Counter = 0;
    for i=Ch
        Counter = Counter + 1;
        Temp1 = Temp([1:ChInfo.nr(i)]+Index(i),:);
        Temp1 = Temp1(:)';
        Data(Counter,[1:length(Temp1)]) = Temp1;
    end
end

Ch = Ch(1);
Slope  = (ChInfo.PhyMax(Ch)-ChInfo.PhyMin(Ch))/(ChInfo.DiMax(Ch)-ChInfo.DiMin(Ch));
Data = (Data-ChInfo.DiMin(Ch))*Slope + ChInfo.PhyMin(Ch);

if or(PlotMode == 3 , PlotMode == 5)
    Data(1,:)=Data(2,:);
end




% --- Executes on button press in CheckBoxAnn.
function CheckBoxAnn_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxAnn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxAnn

if get(hObject,'value')
    set(handles.CheckBoxClaiper,'value',0);
end

if get(hObject,'value')
    set(handles.ListBoxAnnType,'enable','on')
    set(handles.EditAnn,'enable','on')
else
    set(handles.ListBoxAnnType,'enable','off')
    set(handles.EditAnn,'enable','off')
end

% --- Executes on selection change in ListBoxAnnType.
function ListBoxAnnType_Callback(hObject, eventdata, handles)
% hObject    handle to ListBoxAnnType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListBoxAnnType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListBoxAnnType

Temp = get(handles.EditAnn,'string');

if ~isempty(Temp)
    set(handles.EditAnn,'string','');
    Temp = get(hObject,'string');
    Ann=Temp{get(hObject,'value')};
    
    % search for the Ann
    
    for i=1:length(handles.AnnType)
        if strcmp(handles.AnnType{i},Ann);
            Sel = i;
            break
        end
    end
    
    Temp = [];
    Counter = 0;
    for i=[1:length(handles.AnnType)]
        Counter = Counter + 1;
        Temp{Counter} = handles.AnnType{i};
    end
    
    set(handles.ListBoxAnnType,'string',Temp,'value',Sel);
    
end

% --- Executes on selection change in ListBoxAnn.
function ListBoxAnn_Callback(hObject, eventdata, handles)
% hObject    handle to ListBoxAnn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListBoxAnn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListBoxAnn


Sel = get(hObject,'value');

if length(Sel)==1
    Temp = get(handles.PopMenuWindowTime,'string');
    Temp = Temp{get(handles.PopMenuWindowTime,'value')};
    Temp([-3:0]+end)=[];
    WindowTime = str2num(Temp);
    
    
    
    if handles.PlotMode>3
        Temp = (mean(handles.Ann.Time(Sel,1))-WindowTime/4)*handles.SamplingRate;
    else
        
        Temp = (mean(handles.Ann.Time(Sel,1))-WindowTime/4)*handles.SamplingRate;
    end
    
    
    if Temp>get(handles.SliderTime,'max')
        Temp = get(handles.SliderTime,'max');
    end
    if Temp<0
        Temp=0;
    end
    
    set(handles.SliderTime,'value',fix(Temp));
    
    handles=UpdatePlot(handles);
    guidata(hObject,handles);
    
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Temp = get(handles.PopMenuWindowTime,'string');
Temp = Temp{get(handles.PopMenuWindowTime,'value')};
Temp([-3:0]+end)=[];
WindowTime = str2num(Temp);

Sen = get(handles.PopMenuSen,'string');
Sen = Sen{get(handles.PopMenuSen,'value')};
Sen([-2:0]+end)=[];
Sen = str2num(Sen)/2;


if isempty(WindowTime)
    WindowTime = str2num(get(handles.EditCustomSize,'string'));
end


% -------------------------------------------------------

%
Loc=get(handles.axes4,'CurrentPoint');
xLim = get(handles.axes4,'xlim');
yLim = get(handles.axes4,'ylim');

% Check if button down is in histogram of entire R-R window
if (Loc(3)>yLim(1) & Loc(3)<yLim(2) & Loc(1)>xLim(1) & Loc(1)<xLim(2))
    %handles.FlagSaveFile = 1;
    
    
    St = fix(Loc(1)/10)*10;
    
    if get(handles.CheckBoxMinMax,'value')
        [Temp Index]=min(abs(handles.MinMaxValue-St));
        handles.MinMaxValue(Index)=St;
        set(handles.PlotMinMax(Index),'XData',[St St]);
        yLim = get(handles.axes4,'ylim');
        if Index == 1
            set(handles.PlotMinMax(Index+2),'position',[handles.MinMaxValue(1)-200,yLim(2)-diff(yLim)/10],'string',num2str(St));
        else
            set(handles.PlotMinMax(Index+2),'position',[handles.MinMaxValue(2)+20,yLim(2)-diff(yLim)/10],'string',num2str(St));
        end
    end
    
    Index = find(handles.Tac>=St & handles.Tac<(St+10));
    
    if ~isempty(Index)
        Sel = Index(fix(rand*length(Index))+1);
        Shift = fix(handles.IndexTac(Sel)-WindowTime/2*handles.SamplingRate);
        set(handles.SliderTime,'value',Shift);
        
    end
    
    handles = FindHRV_Epoch(handles);
    handles=UpdatePlot(handles);
    guidata(hObject,handles)
end


% if it's not in trace mode
if ~get(handles.CheckBoxTrace,'value')
    
    % ------------------------------
    Loc=get(handles.axes3,'CurrentPoint');
    xLim = get(handles.axes3,'xlim');
    yLim = get(handles.axes3,'ylim');
    
    % Check if button down is in histogram window
    % Click to go to sleep point of sleep axes
    if (Loc(3)>yLim(1) & Loc(3)<yLim(2) & Loc(1)>xLim(1) & Loc(1)<xLim(2))
        Shift = get(handles.SliderTime,'value') + fix((Loc(1)-WindowTime/2)*handles.SamplingRate);
        Temp = get(handles.SliderTime,'max');
        if Shift>Temp
            Shift = Temp;
        end
        
        set(handles.SliderTime,'value',Shift);
        handles=UpdatePlot(handles);
        guidata(hObject,handles)
    end
end



% -----------------------------------------------

Loc=get(handles.axes1,'CurrentPoint');
xLim = get(handles.axes1,'xlim');
yLim = get(handles.axes1,'ylim');

% Check if button down is in histogram window
% Click to go to sleep point of sleep axes
if (Loc(3)>yLim(1) & Loc(3)<yLim(2) & Loc(1)>xLim(1) & Loc(1)<xLim(2))
    
    if get(handles.CheckBoxTrace,'value')
        
        
        Temp = fix(abs(Loc(3)-Sen)/Sen/2);
        
        Shift = get(handles.SliderTime,'value') + fix((Loc(1)-WindowTime/4+WindowTime*Temp)*handles.SamplingRate);
        if Shift<0
            Shift = 0;
        end
        %Shift = get(handles.SliderTime,'value') + WindowTime*Temp*handles.SamplingRate;
        set(handles.SliderTime,'value',Shift);
        set(handles.CheckBoxTrace,'value',0);
        CheckBoxTrace_Callback(handles.CheckBoxTrace, eventdata, handles)
        
    end
    
    %----------------------------------------------
    if get(handles.CheckBoxAnn,'value')
        
        if ~handles.FlagAnn
            Temp = abs(diff(yLim))/100;
            TempX = [Loc(1) Loc(1) Loc(1) Loc(1)];
            TempY = [yLim(1)+Temp yLim(1)+Temp yLim(2)-Temp yLim(2)-Temp];
            
            set(handles.PlotAnn(1),'xData',TempX,'yData',TempY);
            if handles.PlotMode == 1
                set(handles.PlotAnn(2),'xData',TempX,'yData',TempY)
            end
            set(handles.TextAnn,'string','')
            handles.FlagAnn = 1;
        else
            TempX = get(handles.PlotAnn(1),'xData');
            TempX(2:3)=Loc(1);
            set(handles.PlotAnn(1),'xData',TempX);
            if handles.PlotMode == 1
                set(handles.PlotAnn(2),'xData',TempX);
            end
            
            Temp = get(handles.EditAnn,'String');
            
            if isempty(Temp)
                AnnType = get(handles.ListBoxAnnType,'string');
                Ann = AnnType{get(handles.ListBoxAnnType,'value')};
            else
                Ann = Temp;
            end
            
            if isempty(handles.UserAnn)
                Len = 0;
            else
                Len = length(handles.UserAnn.Text);
            end
            handles.UserAnn.Text{Len+1} = Ann;
            handles.UserAnn.Time(Len+1,:) = sort(TempX(1:2)+get(handles.SliderTime,'value')/handles.SamplingRate);
            
            
            Temp = abs(diff(yLim))/10;
            set(handles.TextAnn,'string',Ann,'Position',[min(TempX) yLim(2)-Temp 0]);
            handles.FlagAnn = 0;
            
            % sort Ann with respect to start time
            [Temp Index] = sort(handles.Ann.Time(:,1));
            Counter = 0;
            Temp=[];
            for i=Index'
                Counter = Counter + 1;
                Temp.Text{Counter} = handles.Ann.Text{i};
                Temp.Time(Counter,:) = handles.Ann.Time(i,:);
            end
            handles.Ann = Temp;
            
            %%
            % New annotation has been added and calculate the new HRV 5
            % minutes Index and exclude the new artifact
            %handles.FlagSaveFile = 1;
            handles = FindHRV_Epoch(handles);
            %handles.FlagSaveFile
            
        end
        guidata(hObject,handles)
    end
    
    %-------------------------------------------
    if get(handles.CheckBoxClaiper,'value')
        if ~handles.FlagCliper
            Temp = abs(diff(yLim))/100;
            TempX = [Loc(1) Loc(1) Loc(1) Loc(1)];
            TempY = [yLim(1)+Temp yLim(1)+Temp yLim(2)-Temp yLim(2)-Temp];
            
            set(handles.PlotAnn(1),'xData',TempX,'yData',TempY)
            if handles.PlotMode == 1
                set(handles.PlotAnn(2),'xData',TempX,'yData',TempY)
            end
            set(handles.TextAnn,'string','')
            handles.FlagCliper = 1;
        else
            handles.FlagCliper = 0;
            axes(handles.axes1)
            hold on
            handles.PlotAnn(1) = fill([0 0 0 0], [0 0 0 0] ...
                ,[255 100 50]/255,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
            handles.TextAnn = text(0,0,'','FontWeight','bold','FontSize',12);
            if handles.PlotMode == 1
                axes(handles.axes2)
                hold on
                handles.PlotAnn(2) = fill([0 0 0 0], [0 0 0 0] ...
                    ,[255 100 50]/255,'EdgeColor', [1 1 1],'FaceAlpha',0.25);
            end
        end
        guidata(hObject,handles)
        
    end
    %---------------------------------------------
    
end



% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



Loc=get(handles.axes1,'CurrentPoint');

xLim = get(handles.axes1,'xlim');
yLim = get(handles.axes1,'ylim');

% Check if button down is in histogram window
% Click to go to sleep point of sleep axes
if (Loc(3)>yLim(1) & Loc(3)<yLim(2) & Loc(1)>xLim(1) & Loc(1)<xLim(2))
    
    if get(handles.CheckBoxAnn,'value')
        if handles.FlagAnn
            TempX = get(handles.PlotAnn(1),'xData');
            TempX(2:3)=Loc(1);
            set(handles.PlotAnn(1),'xData',TempX)
            if handles.PlotMode == 1
                set(handles.PlotAnn(2),'xData',TempX)
            end
        end
        guidata(hObject,handles)
    end
    
    if get(handles.CheckBoxClaiper,'value')
        
        if handles.FlagCliper
            TempX = get(handles.PlotAnn(1),'xData');
            TempX(2:3)=Loc(1);
            set(handles.PlotAnn(1),'xData',TempX)
            if handles.PlotMode == 1
                set(handles.PlotAnn(2),'xData',TempX)
            end
            if TempX(1)~=0
                yLim = get(handles.axes1,'yLim');
                Temp = abs(diff(yLim))/5;
                Ann = num2str(fix(abs(diff(TempX(1:2)))*1000));
                set(handles.TextAnn,'string',[Ann],'Position',[min(TempX) yLim(2)-Temp 0]);
                guidata(hObject,handles)
            end
            
        end
    end
    
end




% --- Executes on selection change in PopMenuSen.
function PopMenuSen_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenuSen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenuSen contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenuSen

Sel = get(hObject,'value');
if get(handles.CheckBoxTrace,'value')
    % Trace mode
    handles.SenSel(2) = Sel;
else
    % regular mode
    handles.SenSel(1) = Sel;
end

handles=UpdatePlot(handles);
guidata(hObject,handles);



% --- Executes on button press in RadioButtonGrid.
function RadioButtonGrid_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RadioButtonGrid
handles=UpdatePlot(handles);
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

currentCharacter = get(hObject,'CurrentCharacter')+0;


WindowTime = str2num(get(handles.EditSkipTime,'string'));

Loc=get(handles.axes1,'CurrentPoint');

xLim = get(handles.axes1,'xlim');
yLim = get(handles.axes1,'ylim');


if ~isempty(currentCharacter)
    
    
    switch currentCharacter
        case 32 % space
            if ~get(handles.CheckBoxTrace,'value')
                uicontrol(handles.EditAnn);
                set(handles.CheckBoxAnn,'value',1);
                CheckBoxAnn_Callback(handles.CheckBoxAnn, eventdata, handles);
            end
            
        case 84 % change to trace mode
            
            if get(handles.CheckBoxTrace,'value')
                set(handles.CheckBoxTrace,'value',0);
            else
                set(handles.CheckBoxTrace,'value',1);
            end
            CheckBoxTrace_Callback(handles.CheckBoxTrace, eventdata, handles);
            
            
        case 29
            Temp = get(handles.SliderTime,'value');
            Temp = WindowTime*handles.SamplingRate + Temp;
            Temp1 = get(handles.SliderTime,'max');
            if Temp<Temp1
                set(handles.SliderTime,'value',Temp);
            end
            handles = UpdatePlot(handles);
            guidata(hObject,handles)
            
        case 28
            Temp = get(handles.SliderTime,'value');
            Temp = Temp - WindowTime*handles.SamplingRate ;
            
            if Temp>=0
                set(handles.SliderTime,'value',Temp);
            end
            handles = UpdatePlot(handles);
            guidata(hObject,handles);
            
        case 30
            if (Loc(3)>yLim(1) & Loc(3)<yLim(2) & Loc(1)>xLim(1) & Loc(1)<xLim(2))
                
                Temp = get(handles.PopMenuWindowTime,'Value');
                if Temp<6
                    set(handles.PopMenuWindowTime,'Value',Temp+1);
                    PopMenuWindowTime_Callback(handles.PopMenuWindowTime, eventdata, handles);
                end
            end
        case 31
            if (Loc(3)>yLim(1) & Loc(3)<yLim(2) & Loc(1)>xLim(1) & Loc(1)<xLim(2))
                Temp = get(handles.PopMenuWindowTime,'Value');
                if Temp>1
                    set(handles.PopMenuWindowTime,'Value',Temp-1);
                    PopMenuWindowTime_Callback(handles.PopMenuWindowTime, eventdata, handles);
                end
            end
            
    end
    
end



% --- Executes on button press in PushButtonRemove.
function PushButtonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to PushButtonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Sel = get(handles.ListBoxAnn,'value');

IndexUserAnn = Sel - (length(handles.Ann.Time)-length(handles.UserAnn.Text));

%handles.FlagSaveFile = 1;

if max(Sel) == length(handles.Ann.Text)
    set(handles.ListBoxAnn,'value',Sel-1);
end

Index = [1:length(handles.Ann.Text)];
Index(Sel) = [];

Counter = 0;
Temp=[];
for i=Index
    Counter = Counter + 1;
    Temp.Text{Counter} = handles.Ann.Text{i};
    Temp.Time(Counter,:) = handles.Ann.Time(i,:);
end
handles.Ann = Temp;

Temp=[];
if ~isempty(handles.Ann)
    
    Temp1 = datestr(handles.Ann.Time(:,1)/86400+handles.StartTime,'HH:MM:SS');
    for i=1:length(handles.Ann.Text)
        Temp{i} = [num2str(i) '. ' Temp1(i,:) ...
            '  ' handles.Ann.Text{i}];
    end
    
    
else
    Temp = [];
end
set(handles.ListBoxAnn,'string',Temp);


if IndexUserAnn > 0
    % user ann had been selected and need to removed
    Index = [1:length(handles.UserAnn.Text)];
    Index(IndexUserAnn) = [];
    
    Temp = [];
    Counter = 0;
    
    for i=Index
        Counter = Counter + 1;
        Temp.Text{Counter} = handles.UserAnn.Text{i};
        Temp.Time(Counter,:) = handles.UserAnn.Time(i,:);
    end
    
    handles.UserAnn = Temp;
end



handles=UpdatePlot(handles);

guidata(hObject,handles);


% --- Executes on button press in PushButtonHRV.
function PushButtonHRV_Callback(hObject, eventdata, handles)
% hObject    handle to PushButtonHRV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find the HRV epoch

if 0 % it can be usedc to process everything
    
    FolderIndex = 1:length(handles.FolderNames);
    
    Temp = 'TotalResult.txt';
    
else
    FolderIndex = handles.SelectedFolder;
    
    Temp = [handles.FolderNames{handles.SelectedFolder} '_HRV.txt'];
    
end

[Name Path]=uiputfile(Temp);
FileName = [Path Name];


Fid = fopen(FileName,'w');

Labels = {'ID' 'Epoch no','Start Date/Time','MNN','SDNN','RMSSD','CV','SD1','SD2','SDRatio','LFP','HFP','LHR',...
    'DFA-Alpha1','DFA-Alpha2'};
Temp = [];
for i=1:length(Labels)
    Temp = [Temp Labels{1,i} char(9)];
end
Labels = [Temp char([13 10])];
fwrite (Fid,Labels);
fclose(Fid)


FileIndex = get(handles.ListBoxFileNames,'value');
for CounterFolder = FolderIndex
    
    %set(handles.ListBoxFolderNames,'value',CounterFolder)
    
    
    %handles = ListBoxFolderNames_Callback(handles.ListBoxFolderNames, eventdata, handles);
    
    if get(handles.CheckBoxAllFiles,'value')
        FileIndex = 1:length(handles.FileNames);
    end
        
    
    for CounterFile = FileIndex
        
        set(handles.ListBoxFileNames,'value',CounterFile)
        try
            if get(handles.CheckBoxAllFiles,'value')
                handles = ListBoxFileNames_Callback(handles.ListBoxFileNames, eventdata, handles);
            end
            
            pause(0.1)
            EpochLength = str2num(get(handles.EditHRV_Len,'string')) * handles.SamplingRate;
            
            Sel = get(handles.PopMenuHRV_Ch,'value');
            
            IndexEpoch = FindHRV_Segments(handles.SamplingRate,get(handles.CheckBoxArtifact,'value'),get(handles.CheckBoxMinMax,'value'),handles.R_Index{Sel},...
                handles.Art_Index{Sel},handles.UserAnn,handles.MinMaxValue,EpochLength);
            
            for Counter = 1:length(IndexEpoch)
                % Find the R points
                R_Index = handles.R_Index{Sel}(find((handles.R_Index{Sel} > IndexEpoch(Counter)) & (handles.R_Index{Sel} < (IndexEpoch(Counter)+EpochLength))));
                Time = (R_Index - R_Index(1))/handles.SamplingRate;
                % calculate R-R
                R_R = diff(R_Index)/handles.SamplingRate;
                [SD1C SD2C SD_Ratio CDOWN CUP]=PoinCareFnct(R_R(1:end-1),R_R(2:end),[]);
                [WK1 WK2 JMAX PROB]=FASPER(Time,R_R);
                
                HF_Range = [0.15 0.4];
                LF_Range = [0.04 0.15];
                ULF      = 0.003;
                Total_Range = 0.4;
                
                Index = find(WK1<LF_Range(1));
                WK2(Index) = 0;
                
                Temp = HF_Range;
                Index = find(WK1>Temp(1) & WK1<Temp(2));
                HF_Power = mean(WK2(Index));
                
                Temp = LF_Range;
                Index = find(WK1>Temp(1) & WK1<Temp(2));
                LF_Power = mean(WK2(Index));
                
                SlopeEst = DFA(R_R);
                
                Res=[mean(R_R) std(R_R) sqrt(mean(R_R.^2)) std(R_R)/mean(R_R) SD1C SD2C SD_Ratio LF_Power HF_Power LF_Power/HF_Power SlopeEst];
                
                Text = [handles.FolderNames{handles.SelectedFolder} char(9)...
                    num2str(Counter) char(9) ...
                    datestr(IndexEpoch(Counter)/handles.SamplingRate/86400+handles.StartTime,'mm/dd/yy HH:MM:SS') ...
                    char(9) num2str(Res,['%.4f' char(9)]) char([13 10])];
                
                Fid = fopen(FileName,'a');
                fprintf(Fid,Text);
            end
            
        catch
            [CounterFolder CounterFile]

        end
    end
end

fclose('all');

%%
function EditSkipTime_Callback(hObject, eventdata, handles)
% hObject    handle to EditSkipTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSkipTime as text
%        str2double(get(hObject,'String')) returns contents of EditSkipTime as a double


% Temp = get(handles.PopMenuWindowTime,'string');
% Temp = Temp{get(handles.PopMenuWindowTime,'value')};
% Temp([-3:0]+end)=[];
% WindowTime = str2num(Temp);
%
% if isempty(WindowTime)
%     WindowTime = str2num(get(handles.EditCustomSize,'string'));
%     set(handles.EditCustomSize,'visible','on');
% else
%     set(handles.EditCustomSize,'visible','off');
% end
%
% if WindowTime>str2num(get(handles.EditSkipTime,'string'))
%     set(handles.EditSkipTime,'string',num2str(WindowTime))
% end


% --- Executes on button press in CheckBoxClaiper.
function CheckBoxClaiper_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxClaiper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxClaiper

set(handles.CheckBoxAnn,'value',0);
set(handles.EditAnn,'enable','off');



function EditAnn_Callback(hObject, eventdata, handles)
% hObject    handle to EditAnn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditAnn as text
%        str2double(get(hObject,'String')) returns contents of EditAnn as a double

Key = get(hObject,'string');

if ~isempty(Key)
    Key = Convert2Small(Key);
    
    % search to find the keyword
    Index = [];
    for i=1:length(handles.AnnType)
        Temp = strfind(Convert2Small(handles.AnnType{i}),Key);
        if ~isempty(Temp)
            Index = [Index i];
        end
    end
else
    Index = [1:length(handles.AnnType)];
end

Temp = [];
Counter = 0;
for i=Index
    Counter = Counter + 1;
    Temp{Counter} = handles.AnnType{i};
end

set(handles.ListBoxAnnType,'string',Temp,'value',1);


function Out = Convert2Small(Text)
% convert to small letters
Out = Text;
Temp1 = find(Text>=65 & Text<=90);
Out(Temp1) = Out(Temp1) + 32;


% --- Executes on button press in CheckBoxCh1.
function CheckBoxCh1_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxCh1

% plot Mode
% 1 : Both Ch1 and Ch2
% 2 : Just Ch1
% 3 : Just Ch2
% 4 : Trace plot ch1
% 5 : Trace plot ch2

if get(handles.CheckBoxTrace,'value')
    if get(hObject,'value')
        set(handles.CheckBoxCh2,'value',0);
        handles.PlotMode = 4;
    else
        set(handles.CheckBoxCh2,'value',1);
        handles.PlotMode = 5;
    end
    
    handles.OldPlotMode.Trace = handles.PlotMode;
    
else
    
    
    Sel1 = get(handles.CheckBoxCh1,'value');
    Sel2 = get(handles.CheckBoxCh2,'value');
    
    if Sel1*Sel2
        set(handles.axes2,'visible','on');
        set(handles.axes1,'position',handles.AxesPosition);
        handles.PlotMode = 1;
    else
        axes(handles.axes2)
        cla
        set(handles.axes2,'visible','off');
        Temp = handles.AxesPosition.*[1 1 1 2];
        Temp(2) = Temp(2)-Temp(4)/2;
        set(handles.axes1,'position',Temp);
        handles.PlotMode = 3;
    end
    if Sel1 == 0 & Sel2 ==0
        set(handles.CheckBoxCh2,'value',1);
        handles.PlotMode = 3;
    end
    
    handles.OldPlotMode.Signals = handles.PlotMode;
end

handles=UpdatePlot(handles);

guidata(hObject,handles);


% --- Executes on button press in CheckBoxCh2.
function CheckBoxCh2_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxCh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxCh2


% plot Mode
% 1 : Both Ch1 and Ch2
% 2 : Just Ch1
% 3 : Just Ch2
% 4 : Trace plot ch1
% 5 : Trace plot ch2

if get(handles.CheckBoxTrace,'value')
    if get(hObject,'value')
        set(handles.CheckBoxCh1,'value',0)
        handles.PlotMode = 5;
    else
        set(handles.CheckBoxCh1,'value',1)
        handles.PlotMode = 4;
    end
    
    handles.OldPlotMode.Trace = handles.PlotMode;
else
    
    
    Sel1 = get(handles.CheckBoxCh1,'value');
    Sel2 = get(handles.CheckBoxCh2,'value');
    
    if Sel1*Sel2
        set(handles.axes2,'visible','on');
        set(handles.axes1,'position',handles.AxesPosition);
        handles.PlotMode = 1;
    else
        axes(handles.axes2)
        cla
        set(handles.axes2,'visible','off');
        Temp = handles.AxesPosition.*[1 1 1 2];
        Temp(2) = Temp(2)-Temp(4)/2;
        set(handles.axes1,'position',Temp);
        handles.PlotMode = 2;
    end
    
    if Sel1 == 0 & Sel2 ==0
        set(handles.CheckBoxCh1,'value',1);
        handles.PlotMode = 2;
    end
    
    handles.OldPlotMode.Signals = handles.PlotMode;
end

handles=UpdatePlot(handles);
guidata(hObject,handles);


% --- Executes on button press in CheckBoxTrace.
function CheckBoxTrace_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxTrace


% plot Mode
% 1 : Both Ch1 and Ch2
% 2 : Just Ch1
% 3 : Just Ch2
% 4 : trace plot ch1
% 5 : trace plot ch2

if get(hObject,'value')
    
    
    set(handles.PopMenuSen,'value',handles.SenSel(2));
    set(handles.CheckBoxAnn,'enable','off','value',0)
    set(handles.CheckBoxClaiper,'enable','off','value',0)
    set(handles.ListBoxAnnType,'enable','off')
    set(handles.EditAnn,'enable','off')
    
    set(handles.EditNumTrace,'enable','on');
    
    axes(handles.axes2)
    cla
    set(handles.axes2,'visible','off');
    
    %     axes(handles.axes3)
    %     cla
    %     set(handles.axes3,'visible','off');
    
    Temp = handles.AxesPosition.*[1 1 1 2];
    Temp(2) = Temp(2)-Temp(4)/2;
    
    
    % for removing the axes3 and take the whole screen
    %Temp = handles.AxesPosition.*[1 1 1 2.5];
    %Temp(2) = Temp(2)-Temp(4)/2.5;
    
    
    set(handles.axes1,'position',Temp);
    
    
    if handles.OldPlotMode.Trace ==4
        handles.PlotMode = 4;
        set(handles.CheckBoxCh1,'value',1);
        set(handles.CheckBoxCh2,'value',0);
    else
        handles.PlotMode = 5;
        set(handles.CheckBoxCh1,'value',0);
        set(handles.CheckBoxCh2,'value',1);
    end
else
    set(handles.PopMenuSen,'value',handles.SenSel(1));
    set(handles.CheckBoxAnn,'enable','on')
    set(handles.CheckBoxClaiper,'enable','on')
    
    set(handles.axes1,'position',handles.AxesPosition);
    set(handles.axes2,'visible','on');
    set(handles.axes3,'visible','on');
    
    switch handles.OldPlotMode.Signals
        
        case 1
            set(handles.CheckBoxCh1,'value',1);
            set(handles.CheckBoxCh2,'value',1);
            handles.PlotMode = 1;
        case 2
            set(handles.CheckBoxCh1,'value',1);
            set(handles.CheckBoxCh2,'value',0);
            handles.PlotMode = 2;
        case 3
            set(handles.CheckBoxCh1,'value',0);
            set(handles.CheckBoxCh2,'value',1);
            handles.PlotMode = 3;
    end
    
    set(handles.EditNumTrace,'enable','off');
end

guidata(hObject,handles);

CheckBoxCh1_Callback(handles.CheckBoxCh1, eventdata, handles)



function EditNumTrace_Callback(hObject, eventdata, handles)
% hObject    handle to EditNumTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditNumTrace as text
%        str2double(get(hObject,'String')) returns contents of EditNumTrace as a double



handles=UpdatePlot(handles);


% --- Executes on button press in CheckBoxFilter.
function CheckBoxFilter_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxFilter

handles=UpdatePlot(handles);


% --- Executes on button press in CheckBoxClipping.
function CheckBoxClipping_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxClipping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxClipping


handles=UpdatePlot(handles);


% --- Executes on button press in CheckBoxNotchFilter.
function CheckBoxNotchFilter_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxNotchFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxNotchFilter

handles=UpdatePlot(handles);


% --- Executes on selection change in PopMenuHRV_Ch.
function PopMenuHRV_Ch_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenuHRV_Ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenuHRV_Ch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenuHRV_Ch

handles=UpdatePlot(handles);


% --- Executes on selection change in ListBoxChSelect.
function ListBoxChSelect_Callback(hObject, eventdata, handles)
% hObject    handle to ListBoxChSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListBoxChSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListBoxChSelect

handles=UpdatePlot(handles);





% --- Executes on button press in CheckBoxMinMax.
function CheckBoxMinMax_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxMinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxMinMax

if get(hObject,'value')
    
    yLim=get(handles.axes4,'ylim');
    axes(handles.axes4)
    hold on
    handles.PlotMinMax(1)=plot([1 1]*handles.MinMaxValue(1),yLim,'color','r','linewidth',1.5);
    handles.PlotMinMax(2)=plot([1 1]*handles.MinMaxValue(2),yLim,'color','r','linewidth',1.5);
    handles.PlotMinMax(3) = text(handles.MinMaxValue(1)-200,yLim(2)-diff(yLim)/10,num2str(handles.MinMaxValue(1)),'fontweight','bold');
    handles.PlotMinMax(4) = text(handles.MinMaxValue(2)+20,yLim(2)-diff(yLim)/10,num2str(handles.MinMaxValue(2)),'fontweight','bold');
    
else
    delete(handles.PlotMinMax(1));
    delete(handles.PlotMinMax(2));
    delete(handles.PlotMinMax(3));
    delete(handles.PlotMinMax(4));
    
    
end



handles = FindHRV_Epoch(handles);
handles=UpdatePlot(handles);
guidata(hObject,handles)




function handles = FindHRV_Epoch(handles)

% find the HRV epoch

Index = [];

TempCh = get(handles.ListBoxChSelect,'value');

EpochLength = str2num(get(handles.EditHRV_Len,'string')) * handles.SamplingRate;

Index = FindHRV_Segments(handles.SamplingRate,get(handles.CheckBoxArtifact,'value'),get(handles.CheckBoxMinMax,'value'),handles.R_Index{TempCh(1)},...
    handles.Art_Index{TempCh(1)},handles.UserAnn,handles.MinMaxValue,EpochLength);

handles.Ann.Time = [Index'/handles.SamplingRate Index'/handles.SamplingRate+str2num(get(handles.EditHRV_Len,'string'))];
handles.Ann.Text = [];
if ~isempty(Index)
    for i=1:length(Index)
        handles.Ann.Text{i}='Start of 5 Min HRV';
    end
end


if ~isempty(handles.UserAnn)
    Counter = size(handles.Ann.Time,1);
    handles.Ann.Time = [handles.Ann.Time;handles.UserAnn.Time];
    for i=1:length(handles.UserAnn.Text)
        Counter = Counter + 1;
        handles.Ann.Text{Counter}=handles.UserAnn.Text{i};
    end
end



% [Temp Index] = sort(handles.Ann.Time(:,1));
%
% Temp=[];
% Temp1 = [];
% Counter = 0;
% for i=Index'
%     Counter = Counter + 1;
%     Temp(Counter,:) = handles.Ann.Time(i,:);
%     Temp1{Counter}  = handles.Ann.Text{i};
% end
% handles.Ann.Time = Temp;
% handles.Ann.Text = Temp1;

Temp = [];
if isempty(handles.Ann.Time)
    set(handles.ListBoxAnn,'string','','value',1);
else
    Temp1 = datestr(handles.Ann.Time(:,1)/86400+handles.StartTime,'HH:MM:SS');
    
    for i=1:size(handles.Ann.Time,1)
        Temp{i} = [num2str(i) '. ' Temp1(i,:) ...
            '  ' handles.Ann.Text{i}];
    end
    set(handles.ListBoxAnn,'string',Temp,'value',1);
end

handles=UpdatePlot(handles);


% --- Executes on button press in CheckBoxArtifact.
function CheckBoxArtifact_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxArtifact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxArtifact
%handles.FlagSaveFile = 1;
handles = FindHRV_Epoch(handles);
guidata(hObject,handles);


% --- Executes on button press in CheckBoxLowPass.
function CheckBoxLowPass_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxLowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckBoxLowPass

handles=UpdatePlot(handles);




function Index = FindHRV_Segments(SamplingRate,FlagArtifact,FlagMinMax,R_Index,Art_Index,UserAnn,MinMaxValue,EpochLength)

% find the HRV epoch
CandidateIndex = 1;
Index = [];



while CandidateIndex< (R_Index(end)-EpochLength-1)
    
    % check for 10 sec flat
    Start = CandidateIndex;
    End   = CandidateIndex + EpochLength -1;

    if FlagArtifact
        % check to see there is no artifact
        Temp = find((R_Index(Art_Index)>Start) & (R_Index(Art_Index)<End));
    else
        Temp = [];
    end
    
    if isempty(Temp)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check to see there is no user annotation
        if isempty(UserAnn)
            Temp = [];
        else
            Temp1 = fix(UserAnn.Time(:)*SamplingRate);
            Temp = find(Temp1>Start & (Temp1<End));
        end
        if isempty(Temp)
            % consider Min Max
            % find the detced R
            
            
            Temp = find((R_Index>=Start) & (R_Index<End));
            Tac = diff(R_Index(Temp))/SamplingRate*1000;
            
            if FlagMinMax
                Temp1 = find(Tac<=MinMaxValue(1) | Tac>=MinMaxValue(2));
            else
                Temp1 = [];
            end
            
            if isempty(Temp1)
                if ~isempty(Tac)
                    Index = [Index CandidateIndex];
                end
                CandidateIndex = CandidateIndex + EpochLength;
                
            else
                CandidateIndex =  R_Index(Temp(Temp1(end)+1));
            end
            
        else
            
            CandidateIndex =  R_Index(find(R_Index>Temp1(Temp(end)),1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        if Art_Index(Temp(end)) == length(R_Index)
            CandidateIndex = R_Index(Art_Index(Temp(end)))+1;
        else
            CandidateIndex = R_Index(Art_Index(Temp(end))+1)+1;
        end
    end
    
    
end


% --- Executes on selection change in ListBoxFileNames.
function handles = ListBoxFileNames_Callback(hObject, eventdata, handles)
% hObject    handle to ListBoxFileNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListBoxFileNames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListBoxFileNames


%% this segment will save the annotaions before switch to a new file

% if (handles.FlagSaveFile == 1) | ~(get(handles.CheckBoxArtifact,'value'))
%     %warndlg('Save The Annotation Before Changing the File')
%     MinMaxValue=handles.MinMaxValue;
%     UserAnn = handles.UserAnn;
%     FlagArtifact = get(handles.CheckBoxArtifact,'value');
%     
%     
%     if ispc
%         Path = [handles.Path handles.FolderNames{handles.SelectedFolder} '\'];
%     else
%         Path = [handles.Path handles.FolderNames{handles.SelectedFolder} '/'];
%     end
%     
%     FileName  = [Path handles.FileNames{handles.SelectedFile}([1:end-4]) '_Ann'];
%     
%     save(FileName,'MinMaxValue','UserAnn','FlagArtifact');
%     
%     handles.FlagSaveFile = 0;
%     clear FlagArtifact;
%     handles.FlagSaveFile
%     set(handles.CheckBoxArtifact,'value',1);
% end


%%

FileSel = get(handles.ListBoxFileNames,'value');
handles.SelectedFile = FileSel;

if ispc
    FileName = [handles.Path handles.FolderNames{handles.SelectedFolder} '\' handles.FileNames{FileSel}];
else
    FileName = [handles.Path handles.FolderNames{handles.SelectedFolder} '/' handles.FileNames{FileSel}];
end


[handles.FileInfo handles.ChInfo]=EdfInfo(FileName);

fclose('all');



%% this segment will load the previous processing

% % check to load previous ann, detected R, Artifact and etc
% Temp = FileName;
% Temp([-2:0]+end) = 'mat';
% Temp1 = dir(Temp);
% if isempty(Temp1)
%     % Need processing the data to be able to move forward
%     ProcessEDF_Save(FileName,handles.FileInfo,handles.ChInfo);
%     
% end
% 
% load(Temp);

%%

[R_Index Art_Index SelectedCh SamplingRate]=ProcessEDF(FileName,handles.FileInfo,handles.ChInfo,str2num(get(handles.EditHR_Range,'string')));

handles.Art_Index    = Art_Index;
handles.R_Index      = R_Index;
handles.ECG_Ch       = SelectedCh;
handles.SamplingRate = SamplingRate;

guidata(hObject,handles);


handles.StartTime = datenum([handles.FileInfo.StartDate handles.FileInfo.StartTime],'dd.mm.yyHH.MM.SS');
handles.Fid = fopen(FileName);
guidata(hObject,handles);


% filter design
wo = 60/(handles.SamplingRate/2);  bw = wo/35;
[handles.NotchB,handles.NotchA] = iirnotch(wo,bw); % design the notch filter for the given sampling rate
[handles.FilterB handles.FilterA]=butter(1,4/handles.SamplingRate,'high'); % 2 Hz high pass
[handles.FilterB_LowPass handles.FilterA_LowPass] = butter(2,2/handles.SamplingRate*40,'low'); % 40 Hz Low Pass filter
set(handles.CheckBoxLowPass,'value',0);

Temp = [];
for i=1:length(SelectedCh)
    Temp{i} = handles.ChInfo.Labels(SelectedCh(i),:);
    Temp{i}(Temp{i}==32)=[];
end
set(handles.ListBoxChSelect,'string',Temp);
set(handles.PopMenuHRV_Ch,'string',Temp);

% choose the ch with minumum number ArtiFact
Temp =[];
for i=1:length(R_Index)
    Temp = [Temp ;[length(R_Index{i}) length(Art_Index{i})]];
end

[Temp Sel] = sort(Temp(:,2));
set(handles.ListBoxChSelect,'value',Sel([1 2]));
set(handles.PopMenuHRV_Ch,'value',Sel(1));

% Hist Calculation
Temp = handles.R_Index{Sel(1)};
Index = handles.Art_Index{Sel(1)};


% find the segments with arifact at the begining or end
Tac = (diff(Temp)/handles.SamplingRate)*1000;
IndexTac = handles.R_Index{Sel(1)};
Temp1 = find(diff(Index)>100);

Index = [Index Index-1];
Index (Index>length(Tac))=[];
Index = sort(Index);
Index(Index<1)=[];
Index(diff(Index)==0)=[];

Tac(Index)=[];
IndexTac(Index)=[];
handles.Tac = Tac;
handles.IndexTac = IndexTac;
Tac(Tac>2000)=[];
axes(handles.axes4);

cla
hist(Tac,[200:10:2000]);

set(handles.axes4,'fontweight','bold','xlim',[200 2000]);



%% if there is previous ann

% % Check to see previous Ann
% Temp = FileName;
% Temp([-3:0]+end) = [];
% Temp = [Temp '_Ann.mat'];
% Temp1 = dir(Temp);
% 
% if ~isempty(Temp1)
%     load(Temp);
%     
% 
%     handles.MinMaxValue = MinMaxValue;
%     
%     handles.UserAnn = UserAnn;
%     
%     % check the artifact flag
% 
%      set(handles.CheckBoxArtifact,'value',FlagArtifact)
%     
% else
%     
%    handles.MinMaxValue = fix([mean(Tac)-2.5*std(Tac) mean(Tac)+2.5*std(Tac)]/10)*10;
%    handles.UserAnn     = [];
%    set(handles.CheckBoxArtifact,'value',1);
%     
% end

%%


handles.MinMaxValue = fix([mean(Tac)-2.5*std(Tac) mean(Tac)+2.5*std(Tac)]/10)*10;
handles.UserAnn     = [];
set(handles.CheckBoxArtifact,'value',1);

if get(handles.CheckBoxMinMax,'value')
    
    yLim=get(handles.axes4,'ylim');
    hold on
    handles.PlotMinMax(1)=plot([1 1]*handles.MinMaxValue(1),yLim,'color','r','linewidth',1.5);
    handles.PlotMinMax(2)=plot([1 1]*handles.MinMaxValue(2),yLim,'color','r','linewidth',1.5);
    handles.PlotMinMax(3) = text(handles.MinMaxValue(1)-200,yLim(2)-diff(yLim)/10,num2str(handles.MinMaxValue(1)),'fontweight','bold');
    handles.PlotMinMax(4) = text(handles.MinMaxValue(2)+20,yLim(2)-diff(yLim)/10,num2str(handles.MinMaxValue(2)),'fontweight','bold');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Temp = get(handles.PopMenuWindowTime,'string');
Temp = Temp{get(handles.PopMenuWindowTime,'value')};
Temp([-3:0]+end)=[];
WindowTime = str2num(Temp);

if isempty(WindowTime)
    WindowTime = str2num(get(handles.EditCustomSize,'string'));
end


max = fix(handles.FileInfo.DataRecordDuration*handles.FileInfo.NumberDataRecord-WindowTime)*handles.SamplingRate;
set(handles.SliderTime,'max',max,'value',0);

Temp = [0.2 1]*WindowTime*handles.SamplingRate/max;
set(handles.SliderTime,'SliderStep',Temp)

handles = FindHRV_Epoch(handles);
handles=UpdatePlot(handles);
guidata(hObject,handles)



function EditHRV_Len_Callback(hObject, eventdata, handles)
% hObject    handle to EditHRV_Len (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditHRV_Len as text
%        str2double(get(hObject,'String')) returns contents of EditHRV_Len as a double


handles = FindHRV_Epoch(handles);
handles=UpdatePlot(handles);
guidata(hObject,handles);


