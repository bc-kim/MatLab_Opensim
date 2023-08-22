function varargout = XSens_Bioware_converter(varargin)
% XSENS_BIOWARE_CONVERTER MATLAB code for XSens_Bioware_converter.fig
%      XSENS_BIOWARE_CONVERTER, by itself, creates a new XSENS_BIOWARE_CONVERTER or raises the existing
%      singleton*.
%
%      H = XSENS_BIOWARE_CONVERTER returns the handle to a new XSENS_BIOWARE_CONVERTER or the handle to
%      the existing singleton*.
%
%      XSENS_BIOWARE_CONVERTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in XSENS_BIOWARE_CONVERTER.M with the given input arguments.
%
%      XSENS_BIOWARE_CONVERTER('Property','Value',...) creates a new XSENS_BIOWARE_CONVERTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before XSens_Bioware_converter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to XSens_Bioware_converter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help XSens_Bioware_converter

% Last Modified by GUIDE v2.5 11-Sep-2017 15:36:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @XSens_Bioware_converter_OpeningFcn, ...
                   'gui_OutputFcn',  @XSens_Bioware_converter_OutputFcn, ...
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


% --- Executes just before XSens_Bioware_converter is made visible.
function XSens_Bioware_converter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to XSens_Bioware_converter (see VARARGIN)

% Choose default command line output for XSens_Bioware_converter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes XSens_Bioware_converter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = XSens_Bioware_converter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Upload force plate files
% h=msgbox('Select one C3D file from the force plate');
% waitfor(h);
% 
% [filename,PathName] = uigetfile('*.c3d','Select ONE file from your data collection','MultiSelect','off');
% cd(PathName);
% 
% acq = btkReadAcquisition(filename);
% [analogs, analogsInfo] = btkGetAnalogs(acq);

h=msgbox('Select Kistler data');
waitfor(h);

[FileName,PathName] = uigetfile('*.txt','Select MULTIPLE files from Bioware','MultiSelect','on');
cd(PathName);

for files = 1:length(FileName)

[ForcePlatedata] = Kistler_forceplate_import_v2(FileName{1,files},20,inf);
h = warndlg(FileName{1,files},'Check Kistler file!');
waitfor(h);

% %% Define sample rate
% SR = round((length(ForcePlatedata))/((max(ForcePlatedata(:,1)))-(min(ForcePlatedata(:,1)))));

%% Create Y coordinates for CoP
CoP_Y = zeros(length(ForcePlatedata),1);

%% Rotate forces to match ISB convention
% button = questdlg('Did you correct for ISB convention?','title','Yes','No','default');
% if strcmp(button, 'Yes')
%     Forcedatanew = [ForcePlatedata(:,3),ForcePlatedata(:,4),ForcePlatedata(:,2),...
%        ForcePlatedata(:,6), ForcePlatedata(:,7), ForcePlatedata(:,5),...
%        ForcePlatedata(:,10),CoP_Y, ForcePlatedata(:,9)]; 
   
data_Cop_2 = [ForcePlatedata(:,10),CoP_Y,ForcePlatedata(:,9)];  %CoPx,CoPy,CoPz
data_Cop_1 = [];
data_frs_2 = [ForcePlatedata(:,3),ForcePlatedata(:,4),ForcePlatedata(:,2)]; %Fx, Fy, Fz
data_frs_1 = [];
data_mrs_2 = [ForcePlatedata(:,6), ForcePlatedata(:,7), ForcePlatedata(:,5)];   %Mx, My, Mz
data_mrs_1 = [];
data_frs1 = [data_frs_1 data_Cop_1 data_frs_2 data_Cop_2 data_mrs_1 data_mrs_2];
      
% else
%     Forcedatanew = [ForcePlatedata(:,3),ForcePlatedata(:,4),ForcePlatedata(:,2),...
%        ForcePlatedata(:,6), ForcePlatedata(:,7), ForcePlatedata(:,5),...
%        ForcePlatedata(:,10),ForcePlatedata(:,9)];
% end


%% Save .mot file
% names = ['Fx'   'Fy'    'Fz'    'Mx'    'My'    'Mz'    'Copx'  'Copz'];

% acq_rot = btkNewAcquisition(length(names),length(ForcePlatedata));
% btkSetFrequency(acq_rot, SR);
% 
% for i = 1:header_markers
%     btkSetFrameNumber(acq_rot, frames);
%     btkSetPointLabel(acq_rot, i, names{i});
%     btkSetPoint(acq_rot, i, markersrot.(names{i,1}));
% end
% 
% savemot = strrep(FileName{1,files}, '.c3d', '_rotYZ_dynamic.trc');
% btkWriteAcquisition(acq_rot, savemot);

savemot = strrep(FileName{1,files}, '.txt', '_GRF.mot');

novo_frs={savemot;...
    ['nRows=' num2str(length(ForcePlatedata))];...
    'nColumns=10';...
    'endheader'};
novo_frs(5,(1:10)) = {'time','ground_force_vx','ground_force_vy','ground_force_vz',...
    'ground_force_px','ground_force_py','ground_force_pz',...
    'ground_torque_x','ground_torque_y','ground_torque_z'};
    
data_frs_convert = data_frs1;
for kkk1 = 1:size(data_frs_convert,1)
    for kkk2 = 1:size (data_frs_convert,2)
        if isnan(data_frs_convert(kkk1,kkk2))==1
            data_frs_convert(kkk1,kkk2) = 0;
        end
        novo_frs{kkk1+5,kkk2+1} = num2str(data_frs_convert(kkk1,kkk2));
    end
end
% per1 = 1/SR;
% time2=per1:per1:time1;
time2 = ForcePlatedata(:,1)';
for kk11=1:length(time2)
    novo_frs{kk11+5,1} = num2str(time2(kk11));
end
FILE3=fopen(savemot,'w');
for i=1:size(novo_frs,1)
    for j=1:size(novo_frs,2)
        if j>1 && size(novo_frs{i,j},1)~=0
            fprintf(FILE3,'\t');
        end
        if (i==4 || i==5)&& ( size(novo_frs{i,j},1)==0 ) && j>1
            fprintf(FILE3,'\t');
        end
        if ischar(novo_frs{i,j})
            fprintf(FILE3,novo_frs{i,j});
        end
    end
    fprintf(FILE3,'\n');
end
fclose(FILE3);

clear novo_frs

end

%% Close program
h = warndlg('Data conversion completed!');
waitfor(h);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Select whether to scale from T-pose (extract from a motion file) or from a separate motion file
button = questdlg('Did you collect a separate file for the static pose?','title','Yes','No','default');
if strcmp(button, 'Yes')
    
else
%% Label file
prompt2{1,1} = {'Participant name (e.g. P1)'};
dlg_title2 = 'Input';
num_lines2 = 1;
def2 = {'P1'};
subjectname(1,1) = inputdlg(prompt2{1,1},dlg_title2,num_lines2,def2);

%% Import C3D file
h=msgbox('Select one C3D file to extract the t-pose for static calibration');
waitfor(h);

[filename,PathName] = uigetfile('*.c3d','Select ONE file from your data collection','MultiSelect','off');
cd(PathName);

acq = btkReadAcquisition(filename);
markers = btkGetMarkers(acq);
SR = btkGetPointFrequency(acq);

%% Determine offset between torso and X axis (anterior-posterior)
pRightAcromion = markers.pRightAcromion;
pLeftAcromion = markers.pLeftAcromion;
[Torsowidth] = bodyseg(pRightAcromion(1,1),pRightAcromion(1,2),pLeftAcromion(1,1),pLeftAcromion(1,2),50);

Torso_angle = asind((pRightAcromion(1,1)-pLeftAcromion(1,1))./Torsowidth);  %Calculate the rotation of the torso in relation to the X axis

pIJXY = sqrt((markers.pIJ(1,1).^2)+(markers.pIJ(1,2).^2));
pC7SpinalProcessXY = sqrt((markers.pC7SpinalProcess(1,1).^2)+(markers.pC7SpinalProcess(1,2).^2));

if pIJXY>pC7SpinalProcessXY
    R = roty(Torso_angle);
else
    R = roty(Torso_angle+180);
end

%% Rotate Y and Z
names = fieldnames(markers);

for header_markers = 1:length(names)
    temp = markers.(names{header_markers,1});
    temprot = [temp(1,1),temp(1,3),-temp(1,2)];
    
    for frames = 1
        temprot_roty(frames,:) = R*temprot(frames,:)';        
    end
    
    markersrot.(names{header_markers,1}) = temprot_roty;
end

acq_rot = btkNewAcquisition(length(names),length(temp));
btkSetFrequency(acq_rot, SR);

for i = 1:header_markers
    btkSetFrameNumber(acq_rot, frames);
    btkSetPointLabel(acq_rot, i, names{i});
    btkSetPoint(acq_rot, i, markersrot.(names{i,1}));
end

%% Save TRC file - static posture
% savetrc = strrep(subjectname, '.c3d', '_rotYZ_Static.trc');
savetrc = [subjectname{1,1} '_rotYZ_Static.trc'];
btkWriteAcquisition(acq_rot, savetrc);
clear
end

%% Close program
h = warndlg('Data conversion completed!');
waitfor(h);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Block 1 - CREATE .TRC FILES FROM C3D

%% Import TRC header
h=msgbox('Select multiple C3D data files');
waitfor(h);

[filename,PathName] = uigetfile('*.c3d','Select ONE file from your data collection','MultiSelect','on');
cd(PathName);

for files = 1:length(filename)
    h = warndlg(filename{1,files}); waitfor(h);

acq = btkReadAcquisition(filename{1,files});
markers = btkGetMarkers(acq);
SR = btkGetPointFrequency(acq);

%% Remove static trial from file
names = fieldnames(markers);

for header_markers = 1:length(names)
    temp = markers.(names{header_markers,1});
    temprot = temp(2:length(temp),:);  
    markersnew.(names{header_markers,1}) = temprot;
end

%% Determine if the pelvis has to be aligned with the origin of the GCS (e.g. stationary cycling or treadmill running
button = questdlg('Does the task involve stationary cycling or running?','title','Yes','No','default');
if strcmp(button, 'Yes')
    names = fieldnames(markersnew);
    pSacrum = markersnew.pSacrum;
    pSacrumXY = [mean(pSacrum(:,1)),mean(pSacrum(:,2))];
    
    for header_markers = 1:length(names)
        temp = markersnew.(names{header_markers,1});
        tempaligned = [temp(:,1)-pSacrumXY(1,1),temp(:,2)-pSacrumXY(1,2),temp(:,3)];
        markersnew_aligned.(names{header_markers,1}) = tempaligned;
    end
    
else
    markersnew_aligned = markersnew;
end

%% Determine offset between the pelvis and X axis (anterior-posterior) - align with the X-anterior
pRightASI = markersnew_aligned.pRightASI;
pLeftASI = markersnew_aligned.pLeftASI;
[Pelviswidth] = mean(bodyseg(pRightASI(1,1),pRightASI(1,2),pLeftASI(1,1),pLeftASI(1,2),50));

Pelvis_angle = mean(asind((pRightASI(:,1)-pLeftASI(:,1))./Pelviswidth));  %Calculate the rotation of the pelvis in relation to the X axis

% pRightASIXY = sqrt((markersnew.pRightASI(1,1).^2)+(markersnew.pRightASI(1,2).^2));
% pSacrumXY = sqrt((markersnew.pSacrum(1,1).^2)+(markersnew.pSacrum(1,2).^2));

% if pSacrumXY>pRightASIXY
if Pelvis_angle>0
    R = roty(Pelvis_angle);
else
    R = roty(Pelvis_angle+180);
end

%% Rotate Y and Z
names = fieldnames(markersnew_aligned);

for header_markers = 1:length(names)
    temp = markersnew_aligned.(names{header_markers,1});
    temprot = [temp(:,1),temp(:,3),-temp(:,2)];
    
    for frames = 1:length(temp)
        temprot_roty(frames,:) = R*temprot(frames,:)';        
    end
    
    markersrot.(names{header_markers,1}) = temprot_roty;
end

acq_rot = btkNewAcquisition(length(names),length(temp));
btkSetFrequency(acq_rot, SR);

for i = 1:header_markers
    btkSetFrameNumber(acq_rot, frames);
    btkSetPointLabel(acq_rot, i, names{i});
    btkSetPoint(acq_rot, i, markersrot.(names{i,1}));
end

%% Save TRC file - dynamic
savetrc = strrep(filename{1,files}, '.c3d', '_rotYZ_dynamic.trc');
btkWriteAcquisition(acq_rot, savetrc);

end

% %% Block 2 - CREATE .MOT FILES FROM MVNX
% %% Upload XSENS files (*.mvnx)
% [answerdata,directorydata] = uigetfile('*.mvnx','Choose TWO or more files from XSENS to upload','MultiSelect','on');
% cd(directorydata)
% 
% for files = 1:length(answerdata)
%     h = warndlg(answerdata{1,files}); waitfor(h);
%        
%     tree = load_mvnx(answerdata{1,files});
%     
%     %retrieve the data frames from the subject
%     nSamples = length(tree.subject.frames.frame);  
% 
%     for i = 4:nSamples          
%         %Joint Angles
% %        Coronal, transversal, sagital
% %        Coronal = positive for adduction
% %        Transversal = positive for internal rotation
% %        Sagital = Flexion projected angle
%         
% %         UpperSacrum(i,1) = tree.subject.frames.frame(i).jointAngle(1,(1*3)-2);
% %         UpperSacrum(i,2) = tree.subject.frames.frame(i).jointAngle(1,(1*3)-1);   %Id 1 * 3 axes + 1 rotation (Y after X)
% %         UpperSacrum(i,3) = tree.subject.frames.frame(i).jointAngle(1,(1*3));  
% 
%         R_Hip(i,1) = tree.subject.frames.frame(i).jointAngle(1,(15*3)-2);
%         R_Hip(i,2) = tree.subject.frames.frame(i).jointAngle(1,(15*3)-1);   %Id 15 * 3 axes + 1 rotation (Y after X)
%         R_Hip(i,3) = tree.subject.frames.frame(i).jointAngle(1,(15*3));     
%         
%         R_Knee(i,1) = tree.subject.frames.frame(i).jointAngle(1,(16*3)-2);
%         R_Knee(i,2) = tree.subject.frames.frame(i).jointAngle(1,(16*3)-1);   %Id 16 * 3 axes + 1 rotation (Y after X)
%         R_Knee(i,3) = tree.subject.frames.frame(i).jointAngle(1,(16*3)); 
%         
%         R_Ankle(i,1) = tree.subject.frames.frame(i).jointAngle(1,(17*3)-2);
%         R_Ankle(i,2) = tree.subject.frames.frame(i).jointAngle(1,(17*3)-1);   %Id 17 * 3 axes + 1 rotation (Y after X)
%         R_Ankle(i,3) = tree.subject.frames.frame(i).jointAngle(1,(17*3)); 
%         
%         L_Hip(i,1) = tree.subject.frames.frame(i).jointAngle(1,(19*3)-2);
%         L_Hip(i,2) = tree.subject.frames.frame(i).jointAngle(1,(19*3)-1);   %Id 19 * 3 axes + 1 rotation (Y after X)
%         L_Hip(i,3) = tree.subject.frames.frame(i).jointAngle(1,(19*3));     
%         
%         L_Knee(i,1) = tree.subject.frames.frame(i).jointAngle(1,(20*3)-2);
%         L_Knee(i,2) = tree.subject.frames.frame(i).jointAngle(1,(20*3)-1);   %Id 20 * 3 axes + 1 rotation (Y after X)
%         L_Knee(i,3) = tree.subject.frames.frame(i).jointAngle(1,(20*3)); 
%         
%         L_Ankle(i,1) = tree.subject.frames.frame(i).jointAngle(1,(21*3)-2);
%         L_Ankle(i,2) = tree.subject.frames.frame(i).jointAngle(1,(21*3)-1);   %Id 21 * 3 axes + 1 rotation (Y after X)
%         L_Ankle(i,3) = tree.subject.frames.frame(i).jointAngle(1,(21*3)); 
%     end
%     
% %% Save .mat data
% clear h hObject eventdata handles
% filenameOUT = strrep(answerdata{1,files}, '.mvnx', '.mat');
% save(filenameOUT);
% 
% %% Create .mot file
% SR = tree.subject.frameRate;
% Time = (1/SR):(1/SR):((length(L_Ankle))/SR);
% 
% allAnglesMatrix = [Time', R_Hip, R_Knee(:,3), R_Ankle(:,3), L_Hip, L_Knee(:,3), L_Ankle(:,3)];
% [numRows,numCols] = size(allAnglesMatrix);
% ncols = num2str(numCols);
% nrows = num2str(numRows);
% % Create base joint angle matrix
% % splitfilename = strsplit(filename,'.');
% % splitfilename = splitfilename{1};
% filenameOUT = strrep(answerdata{1,files}, '.mvnx', '_Xsens_LowerLimbAngles.mot');
% fileID = fopen(filenameOUT,'w');
% 
% fprintf(fileID,('Coordinates'),'\n');
% fprintf(fileID,strcat(['nnRows=',nrows,'\n\n']));
% fprintf(fileID,strcat(['nnColumns=',ncols,'\n\n']));
% fprintf(fileID,(''),'\n');
% fprintf(fileID,('Units are S.I. units (second, meters, Newtons, ...)'),'\n');
% fprintf(fileID,('Angles are in degrees'),'\n');
% fprintf(fileID,(''),'\n');
% fprintf(fileID,('endheader'),'\n');
% 
% % fprintf(fileID,strcat(['first trial\nnRows=',nrows,'\nnColumns=',ncols,'\n\n']));
% % fprintf(fileID,'# SIMM Motion File Header:\n');
% % nrange = num2str([min(min(allAnglesMatrix)) max(max(allAnglesMatrix))]); 
% % otherdata = num2str(1); %always 1 since time vec exists
% % % jointNameStrings = ['rHip'  'rKnee'  'rAnkleL_Hip'    'lHip' 'lKnee' 'lAnkle'];
% % 
% % fprintf(fileID,strcat(['name XSensLowerLimbAngles\ndatacolumns ',ncols,'\ndatarows ', nrows,'\notherdata ',otherdata,'\nrange ',nrange,'\nendheader\n']));
% % % fprintf(fileID,strcat(['time ',jointNameStrings,'\n']));
% % fprintf(fileID,strcat(['time ','rHip','rKnee','rAnkleL_Hip','lHip','lKnee','lAnkle','\n']));
% % % fprintf(fileID,'%5.3f\t%4.3f\n',Mat);
% fprintf(fileID,'%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t\n',allAnglesMatrix);
% fprintf(fileID,'\n');
% fclose(fileID);
% 
% end

%% Close program
h = warndlg('Data conversion completed!');
waitfor(h);