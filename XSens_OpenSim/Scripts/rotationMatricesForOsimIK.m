%rotationMatricesForOsimIK.m
%Contributors: Johanna O'Day
%Get limb orientations from MVNX file and create file of rotation matrices
%that you can feed into OpenSim IK

%This file takes MVNX file and allows you to visualize the sensor orientations 
%of the thigh and shank sensors with respect to the calculated knee flexion
%angle that Xsens gives you. 
%You then need to run inverse kinematics with these rotation matrices to back out
%your own knee flexion angle from this rotation matrix data.
%%
close all
clear all
%First : Run main_mvnx file
main_mvnx_2

Time = zeros(nSamples,1);
for jj=[1:nSamples]
    Time(jj)=tree.subject.frames.frame(jj).time;
end
Time=Time./1000;
%%
%For automated answers to make this easier for BIOE485 teaching Team
%Comment this out if you are not the BIOE485 teaching team

filenamesplit = strsplit(filename,'bpm');
BPM = strsplit(filenamesplit{1},'_');
BPM = BPM{2};

if str2num(BPM) == 40
    % cropping Xsens data to video for 40 BPM
    Xsens_straight = 20.85; %sec
elseif str2num(BPM) == 80
    %for 80 BPM
    Xsens_straight = 11.97;
elseif str2num(BPM) == 120 
    %for 120 BPM
    Xsens_straight = 9.84; %sec
end

%%
%Crop Xsens data to necessary time
prompt = {'Enter start time:','Enter experiment duration:'};
            dlg_title = 'Cropping Xsens data';
            num_lines = 1;
            defaultans = {num2str(Xsens_straight), num2str(10)}; %default experiment duration is 10 s
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            Xsens_straight = str2num(answer{1,1});
            exp_duration = str2num(answer{2,1});
            
%Determine what index Xsens is at when video starts            
TimeStartIndex = find(Time<=Xsens_straight);
TimeStartIndex = TimeStartIndex(end)
%Determine what index Xsens is at when video ends 
TimeEndIndex = find(Time>=exp_duration+Xsens_straight,1) 
newTime = Time(TimeStartIndex:TimeEndIndex);
%%
%Pull out the separate euler parameters w,x,y,z, from (W,[x,y,z]) for sensor on upper leg and
%lower leg
for sampleIndex = 4:nSamples
    rightUpperLegOrientationW(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(45);
    rightUpperLegOrientationX(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(46);
    rightUpperLegOrientationY(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(47);
    rightUpperLegOrientationZ(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(48);
    rightLowerLegOrientationW(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(49);
    rightLowerLegOrientationX(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(50);
    rightLowerLegOrientationY(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(51);
    rightLowerLegOrientationZ(:,sampleIndex) = tree.subject.frames.frame(sampleIndex).sensorOrientation(52);
    
end
%%
%Make a matrix of the values  [w,x,y,z]
%General case make matrix of upper orientations 
% % rightUpperLegOrientations = [-1*rightUpperLegOrientationW', rightUpperLegOrientationX', rightUpperLegOrientationY', rightUpperLegOrientationZ'];
% % rightLowerLegOrientations = [-1*rightLowerLegOrientationW', rightLowerLegOrientationX', rightLowerLegOrientationY', rightLowerLegOrientationZ'];

%Create rotation matrices from quaternion data (Using function from
%Aerospace toolbox quat2dcm)
%Careful: quat2dcm needs -1*W in qutaernion to work properly???

% rightUpperLegRotMat = quat2dcm(rightUpperLegOrientations);
% rightLowerLegRotMat =  quat2dcm(rightLowerLegOrientations);
%%
%Convert to rotation matrix by equation for rotation matrix transformation
%from quaternions
rightUpperLegOrientations = [rightUpperLegOrientationW', rightUpperLegOrientationX', rightUpperLegOrientationY', rightUpperLegOrientationZ'];
rightLowerLegOrientations = [rightLowerLegOrientationW', rightLowerLegOrientationX', rightLowerLegOrientationY', rightLowerLegOrientationZ'];

%crop orientations to video time
rightUpperLegOrientations = rightUpperLegOrientations(TimeStartIndex:TimeEndIndex,:);
rightLowerLegOrientations = rightLowerLegOrientations(TimeStartIndex:TimeEndIndex,:);

%Index the matrices 

for ii = 1:length(rightUpperLegOrientations)
q = rightUpperLegOrientations(ii,:);
s = norm(q)^-2;
qr=q(1);
qi=q(2);
qj=q(3);
qk=q(4);

rightUpperLegRotMat(:,:,ii) = [ 1-2*s*(qj^2 + qk^2), 2*s*(qi*qj-qk*qr), 2*s*(qi*qk+qj*qr);
    2*s*(qi*qj + qk*qr), 1-2*s*(qi^2 + qk^2), 2*s*(qj*qk - qi*qr);
    2*s*(qi*qk - qj*qr), 2*s*(qj*qk + qi*qr), 1-2*s*(qi^2 + qj^2)];
end

for ii = 1:length(rightLowerLegOrientations)
q = rightLowerLegOrientations(ii,:);
s = norm(q)^-2;
qr=q(1);
qi=q(2);
qj=q(3);
qk=q(4);

rightLowerLegRotMat(:,:,ii) = [ 1-2*s*(qj^2 + qk^2), 2*s*(qi*qj-qk*qr), 2*s*(qi*qk+qj*qr);
    2*s*(qi*qj + qk*qr), 1-2*s*(qi^2 + qk^2), 2*s*(qj*qk - qi*qr);
    2*s*(qi*qk - qj*qr), 2*s*(qj*qk + qi*qr), 1-2*s*(qi^2 + qj^2)];
end

%%
%Save these matrices in file format to read into our inverse kinematics
%script, futureOrientationInverseKinematics.cpp
%Format of text file that is input to futureOrientationInverseKinematics.cpp  
%PacketCounter-SampleTimeFine	Mat[1][1]	Mat[2][1]	Mat[3][1]	Mat[1][2]	Mat[2][2]	Mat[3][2]	Mat[1][3]	Mat[2][3]	Mat[3][3]

%For some reason the Xsens samples don't have any recordings until sample #4
%so start is set at 4
sampleStart = 4;
rightLowerLegRotMat = rightLowerLegRotMat(:,:,sampleStart:end);
rightUpperLegRotMat = rightUpperLegRotMat(:,:,sampleStart:end);

finalRightLowerMat= zeros(length(rightLowerLegRotMat),9);
finalRightUpperMat= zeros(length(rightUpperLegRotMat),9);

for numsamples = 1:length(rightLowerLegRotMat)
    finalRightLowerMat(numsamples,:) = [rightLowerLegRotMat(1,1,numsamples), rightLowerLegRotMat(2,1,numsamples), rightLowerLegRotMat(3,1,numsamples), rightLowerLegRotMat(1,2,numsamples), rightLowerLegRotMat(2,2,numsamples), rightLowerLegRotMat(3,2,numsamples), rightLowerLegRotMat(1,3,numsamples), rightLowerLegRotMat(2,3,numsamples), rightLowerLegRotMat(3,3,numsamples)];
    finalRightUpperMat(numsamples,:) = [rightUpperLegRotMat(1,1,numsamples), rightUpperLegRotMat(2,1,numsamples), rightUpperLegRotMat(3,1,numsamples), rightUpperLegRotMat(1,2,numsamples), rightUpperLegRotMat(2,2,numsamples), rightUpperLegRotMat(3,2,numsamples), rightUpperLegRotMat(1,3,numsamples), rightUpperLegRotMat(2,3,numsamples), rightUpperLegRotMat(3,3,numsamples)];
end

%Packetcounter is filler
packetcounter = zeros(length(rightLowerLegRotMat),1); 
finalRightLowerMat = [packetcounter, finalRightLowerMat];
finalRightUpperMat = [packetcounter, finalRightUpperMat];

% dlmwrite(strcat(filename,'_right_shank_matrix_1.txt'),finalRightLowerMat,'\t');
% dlmwrite(strcat(filename,'_right_thigh_matrix_1.txt'),finalRightUpperMat,'\t');
%%
%For ease, make filenames specific to Melissa's experiment files for BIOE485 teaching
%team, because these are the filenames already in the .cpp file that you're going to next
%If you are not the BIOE485 teaching team then you don't need this
%code
newfilename = BPM;

%%
savepath = pwd;
newShankFileName = fullfile(savepath,strcat(newfilename,'_right_shank_matrix.txt'));
newThighFileName = fullfile(savepath,strcat(newfilename,'_right_thigh_matrix.txt'));

%Throw in file headers to be compatible with XSens Headers
fidShank = fopen(newShankFileName,'wt');
fprintf(fidShank, '// Start Time: Unknown','\n');
fprintf(fidShank,'\n');
fprintf(fidShank,'// Update Rate: 100.0Hz','\n');
fprintf(fidShank,'\n');
fprintf(fidShank,'// Filter Profile: human (46.1)','\n');
fprintf(fidShank,'\n');
fprintf(fidShank,'// Option Flags: AHS Disabled ICC Disabled','\n');
fprintf(fidShank,'\n');
fprintf(fidShank,'// Firmware Version: 4.0.2','\n');
fprintf(fidShank,'\n');
fprintf(fidShank,'PacketCounter	SampleTimeFine	Mat[1][1]	Mat[2][1]	Mat[3][1]	Mat[1][2]	Mat[2][2]	Mat[3][2]	Mat[1][3]	Mat[2][3]	Mat[3][3]','\n');  % header
fprintf(fidShank,'\n');
fclose(fidShank);
dlmwrite(newShankFileName,finalRightLowerMat,'delimiter','\t','-append');

fidThigh = fopen(newThighFileName,'wt');
fprintf(fidThigh, '// Start Time: Unknown','\n');
fprintf(fidThigh,'\n');
fprintf(fidThigh,'// Update Rate: 100.0Hz','\n');
fprintf(fidThigh,'\n');
fprintf(fidThigh,'// Filter Profile: human (46.1)','\n');
fprintf(fidThigh,'\n');
fprintf(fidThigh,'// Option Flags: AHS Disabled ICC Disabled','\n');
fprintf(fidThigh,'\n');
fprintf(fidThigh,'// Firmware Version: 4.0.2','\n');
fprintf(fidThigh,'\n');
fprintf(fidThigh,'PacketCounter	SampleTimeFine	Mat[1][1]	Mat[2][1]	Mat[3][1]	Mat[1][2]	Mat[2][2]	Mat[3][2]	Mat[1][3]	Mat[2][3]	Mat[3][3]','\n');  % header
fclose(fidThigh);
dlmwrite(newThighFileName,finalRightUpperMat,'delimiter','\t','-append');

% dlmwrite(newShankFileName,finalRightLowerMat,'\t');
% dlmwrite(strcat(filename,'back_right_thigh_matrix_1.txt'),finalRightUpperMat,'\t');

%%
%TODO:
%Generalize this code to work for any pairs of sensors