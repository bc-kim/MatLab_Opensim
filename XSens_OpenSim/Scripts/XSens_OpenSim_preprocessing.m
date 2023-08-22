function XSens_OpenSim_preprocessing

% XSens pre-processing for OpenSim
% by Bini
% date August 2017
% -------------------------------------------------------------------------
% Uses .c3d files exported from MVN Studio.
% Extracts one static trial from the data collection to scale the model.
% Rotates coordinates to match with ISB and OpenSim convention.
% Uses the following supporting scripts:
% bodyseg: computes segment length
% btkReadAcquisition
% btkGetMarkers
% =========================================================================

%% BLOCK 1- Extract N-pose and save as static calibration

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

%% BLOCK 2- Import data collection files and rotate axes

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

%% Close program
h = warndlg('Data conversion completed!');
waitfor(h);

end