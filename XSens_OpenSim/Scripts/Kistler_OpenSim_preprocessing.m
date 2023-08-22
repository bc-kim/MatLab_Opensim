function Kistler_OpenSim_preprocessing

% Kistler pre-processing for OpenSim
% by Bini
% date September 2017
% -------------------------------------------------------------------------
% Uses .txt files exported from Bioware.
% Extracts one static trial from the data collection to scale the model.
% Rotates coordinates to match with ISB and OpenSim convention.
% Uses the following supporting scripts:
% Kistler_forceplate_import_v2.m: imports .txt files exported from Bioware
% =========================================================================

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

end