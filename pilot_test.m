
addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

%% Flanker first pilot: Aug 16
datapath = '/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-24531/ses-01/meg/';

bvfiles = {'24531_Flanker3_redone_2022_Aug_16_1521.csv';'24531_Flanker3_redone_2022_Aug_16_1529.csv'};
processingfolder = '/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-24531/ses-01/';

fids_name = 'sub-24531_fiducials.tag';
% 

filenames_old = cell(1,2);
for ii = 1:2
    filenames_old{ii} = ['24531_CognitiveControl_20220816_00',num2str(ii+4),'.ds'];
end
cd(datapath)

filenames = {'sub-24531_ses-01_task-flanker_run-1_meg.ds';'sub-24531_ses-01_task-flanker_run-2_meg.ds'};

% Flanker_process(filenames,datapath,processingfolder,bvfiles,fids_name)

%% Flanker second pilot: Sep 27

datapath = '/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-24531/ses-02/meg/';

bvfiles = {'24531_Flanker3_redone_2022_Sep_27_1719.csv';'24531_Flanker3_redone_2022_Sep_27_1728.csv'};
processingfolder = '/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-24531/ses-02/';
% fids_name = 'sub-24531_fiducials_runs57.tag';
fids_name = 'sub-24531_fiducials.tag';

filenames_old = cell(1,2);
for ii = 1:2
    filenames_old{ii} = ['24531_CognitiveControl_20220927_00',num2str(ii+5),'.ds'];
end
filenames = {'sub-24531_ses-02_task-flanker_run-1_meg.ds';'sub-24531_ses-02_task-flanker_run-2_meg.ds'};

cd(datapath)


% Flanker_process(filenames,datapath,processingfolder,bvfiles,fids_name)

%% AX-CPT first pilot: Aug 16
datapath = '/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-24531/ses-01/meg/';

% bvfiles = {'24531_Flanker3_redone_2022_Aug_16_1521.csv';'24531_Flanker3_redone_2022_Aug_16_1529.csv'};
processingfolder = '/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-24531/ses-01/';

fids_name = 'sub-24531_fiducials.tag';
% 

filenames_old = cell(1,4);
filenames = cell(1,4);
for ii = 1:4
    filenames_old{ii} = ['24531_CognitiveControl_20220816_00',num2str(ii),'.ds'];
    filenames{ii} = ['sub-24531_ses-01_task-axcpt_run-',num2str(ii),'_meg.ds'];
end
cd(datapath)


% AXCPT_process(filenames,datapath,processingfolder,bvfiles,fids_name)

%% AX-CPT second pilot: Sep 27

datapath = '/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-24531/ses-02/meg/';

% bvfiles = {'24531_Flanker3_redone_2022_Sep_27_1719.csv';'24531_Flanker3_redone_2022_Sep_27_1728.csv'};
processingfolder = '/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-24531/ses-02/';
% fids_name = 'sub-24531_fiducials_runs57.tag';
fids_name = 'sub-24531_fiducials.tag';

filenames_old = cell(1,3);
filenames = cell(1,3);

for ii = 1:3
    filenames_old{ii} = ['24531_CognitiveControl_20220927_00',num2str(ii+2),'.ds'];
    filenames{ii} = ['sub-24531_ses-02_task-axcpt_run-',num2str(ii),'_meg.ds'];
end

cd(datapath)


% Flanker_process(filenames,datapath,processingfolder,bvfiles,fids_name)
%%
task = 'axcpt';
for ii =1:3
    sdan = '24531';
    cfg=[];

    cfg.method                  = 'copy';
    % cfg.proc                    = 'raw';
    cfg.dataset                 = filenames_old{ii};
    cfg.bidsroot                = '/data/liuzzil2/MEG_AXCPT_Flanker/data/'; %top level directory for the BIDS output
    cfg.sub                     = sdan; %subject name
    cfg.run                     = ii; %optional
    cfg.ses                     = '02'; %optional
    cfg.datatype                = 'meg'; %can be any of 'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'
    % cfg.participants.age        = meginfo.Age(iiS);
    % cfg.participants.sex        = meginfo.Sex(iiS); %'m' or 'f'
    % cfg.participants.group      = meginfo.Group(iiS);

    cfg.outputfile              = filenames{ii};
    
    if strcmp(task,'flanker')
    cfg.task                    = 'flanker'; %task name is required for functional data
    cfg.TaskName                = cfg.task;
    cfg.TaskDescription         = 'Flanker task, 4 blocks of 48 trials each, congruent/incongruent, left and right stimuli with equal likelihood.';
    else
         cfg.task                    = 'axcpt'; %task name is required for functional data
    cfg.TaskName                = cfg.task;
    cfg.TaskDescription         = 'AX-CPT task, 4 blocks of 80 trials each.';
    end

    cd(cfg.dataset)
    fileID = fopen([cfg.dataset(1:end-2),'infods'],'r');
    TaskDate = [];
    while isempty(TaskDate)
        tline = fscanf(fileID,'%s',1);
    %     tline = fgetl(fileID);
        if contains(tline,'DATASET_COLLECTIONDATETIME')
            tline = fscanf(fileID,'%s',1);

            ind20 = strfind(tline,'20'); % find start of date, i.e. 2019 or 2020
            TaskDate = tline(ind20(1)+[0:13]);
        end
    end
    fclose(fileID);


    cfg.scans.acq_time          =  [TaskDate(1:4),'-',TaskDate(5:6),'-',TaskDate(7:8),...
        'T',TaskDate(9:10),':',TaskDate(11:12),':',TaskDate(13:14)];
    cd ..
    disp(cfg.scans.acq_time)

    data2bids(cfg);
end



