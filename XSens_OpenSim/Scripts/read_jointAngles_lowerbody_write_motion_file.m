%read_jointAngles_lowerbody_write_motion_file.m
%Contributors: Johanna O'Day, Melissa Boswell
% This script reads an .mvnx file collected with Xsens Awinda IMU sensors
% in MVN Studio and takes creates a table of all the joints and joint angle
% estimations at every time point (jointTable)
% These estimates come from Xsens proprietary algorithms and biomechanical model.

%It then outputs selected joint angles to a matrix with time as first row 
%(allAnglesMatrix)
%Finally, it creates a motion file named ?mvnxfilename_Xsens_jointangle_q.mot? 
%of joint angles to input in Open Sim from Xsens IMU Data
%   This file has the required OpenSim Headings 
% and has the following parameters: pelvis_tilt, pelvis_tx, pelvis_ty, hip_flexion_r, knee_angle_r, ankle_angle_r

%Future improvements could generalize these parameters so that you can
%specify which joint angles you want to get from the Xsens file.

%%

% This runs main_mvnx.m (this is code modified from Xsens Toolkit (Xsens North America Inc.)
% main_mvnx.m loads two figures of the ?first segment position? - the pelvis rotation in 3D and the pelvis 3D translation.
% It uses the load_mvnx function (Xsens Toolkit) and creates a data tree
% 'tree'

% main_mvnx
main_mvnx_2

%Create joint table of all joints and joint angles from IMU data

[~,nJoints] = size(tree.subject.joints.joint);
jointNames = cell(3*nJoints+1,1);

for hh = 1:nJoints
    jointName = tree.subject.joints.joint(hh).label;
    %repeat jointName 3 times
    jointNames{hh*3-2,1} = jointName;
    jointNames{hh*3-1,1} = jointName;
    jointNames{hh*3,1} = jointName;
    %add time if you're at the last joint
    if hh == nJoints
        jointNames{hh*3+1,1} = cellstr('time');
    end
end

%%

ColIndx = repmat([1,2,3]',nJoints,1); %These correspond to X,Y,Z angles (ColIndx = 1,2,3) 
ColIndx = [ColIndx; 0]; %add a zero for time

Angles = zeros(length(jointNames),nSamples);  %Angles in degrees
jointTable = table(jointNames,ColIndx,Angles);

for ii = 1:length(jointTable.jointNames)-1
    temp = zeros(1,nSamples);
    for sampleIndex = 4:nSamples
        temp(:,sampleIndex)= tree.subject.frames.frame(sampleIndex).jointAngle(ii); %joint angles
    end
    jointTable.Angles(ii,:) = temp;
end
%%
%Create time vector and add to table     
Time = zeros(nSamples,1);
for jj=[1:nSamples]
    Time(jj)=tree.subject.frames.frame(jj).time;
end
Time=Time./1000;

jointTable.Angles(length(jointNames),:) = Time;

%This gives you a joint table with each joint segment (as labeled in MVN
%file and then gives you a row vector of each X,Y,Z angles for each joint over time)
% Now you can pull out vectors that you want


%%
% Pull out the angles we want for this model
% pelvis_tilt	= Z (+/-)
% pelvis_tx	= x (position 1 x)
% pelvis_ty	= y  (position 2 y)
% hip_flexion_r	= Z (+/-)
% knee_angle_r = only Z	
% ankle_angle_r = Z

%For some reason the Z-axis angles are opposite in OpenSim and Xsens, so we
%multiply by a negative 1 to invert them so that we get correct motion in
%OpenSim
negZ = -1;

%OpenSim names
rows = strcmp(jointTable.jointNames,'jL5S1') & jointTable.ColIndx==3;
vars = {'Angles'};
pelvis_tilt = negZ*jointTable{rows,vars};

%Just find these the hard way for now, not in table? 
%TODO: Add to table later
for sampleIndex = 4:nSamples
    pelvis_tx(:,sampleIndex)= tree.subject.frames.frame(sampleIndex).position(1); %x
    pelvis_ty(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).position(2); %y
end

rows = strcmp(jointTable.jointNames,'jRightHip') & jointTable.ColIndx==3;
vars = {'Angles'};
hip_flexion_r = negZ*jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightKnee') & jointTable.ColIndx==3;
vars = {'Angles'};
knee_angle_r = negZ*jointTable{rows,vars};

rows = strcmp(jointTable.jointNames,'jRightAnkle') & jointTable.ColIndx==3;
vars = {'Angles'};
ankle_angle_r = negZ*jointTable{rows,vars};

%%
%Now output all to matrix

allAnglesMatrix = [Time'; pelvis_tilt; pelvis_tx;  pelvis_ty; hip_flexion_r; knee_angle_r; ankle_angle_r];

%%
% Creates motion file of joint angles 
% to input in Open Sim from XSens IMU Data

[numRows,numCols] = size(allAnglesMatrix');
ncols = num2str(numCols);
nrows = num2str(numRows);
% Create base joint angle matrix
splitfilename = strsplit(filename,'.');
splitfilename = splitfilename{1};
fileID = fopen(strcat(splitfilename,'_Xsens_jointangle_q','.mot'),'w');
fprintf(fileID,strcat(['first trial\nnRows=',nrows,'\nnColumns=',ncols,'\n\n']));
fprintf(fileID,'# SIMM Motion File Header:\n');
nrange = num2str([min(min(allAnglesMatrix)) max(max(allAnglesMatrix))]); 
otherdata = num2str(1); %always 1 since time vec exists
jointNameStrings = ['pelvis_tilt pelvis_tx pelvis_ty	hip_flexion_r knee_angle_r ankle_angle_r'];

fprintf(fileID,strcat(['name Knee Extension\ndatacolumns ',ncols,'\ndatarows ', nrows,'\notherdata ',otherdata,'\nrange ',nrange,'\nendheader\n']));
fprintf(fileID,strcat(['time ',jointNameStrings,'\n']));
% fprintf(fileID,'%5.3f\t%4.3f\n',Mat);
fprintf(fileID,'%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t\n',allAnglesMatrix);
fprintf(fileID,'\n');
fclose(fileID);
