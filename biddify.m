addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath ~/matlab_utils/
%%
sdan = '24663';
datapath = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sdan,'/meg/'];


task = 'axcpt'; 
% recs = 7:9;
% d = dir(datapath);
% filenames_old = cell(1,length(recs));
% filenames = cell(1,length(recs));
% 
% 
% for ii = 1:length(recs)
%     filenames_old{ii} = d(ii+recs(1)+1).name;
%     filenames{ii} = ['sub-',sdan,'_task-',task,'_run-',num2str(ii),'_meg.ds'];
% end


recs = 1:2;

d = dir(datapath);
filenames_old = cell(1,length(recs));
filenames = cell(1,length(recs));
r = 0;
for ii = 6:7%length(d)
%     if contains( d(ii).name ,'Cognitive')
        r = r+1;
        filenames_old{r} = d(ii).name;
        filenames{r} = ['sub-',sdan,'_task-',task,'_run-',num2str(r),'_meg.ds'];
%         filenames{r} = [d(ii).name(1:end-5), '.ds'];
%     end
end

cd(datapath)



%%

for ii =1:length(filenames_old)
    
    cfg=[];

    
    % cfg.proc                    = 'raw';
    cfg.dataset                 = filenames_old{ii};
    cfg.bidsroot                = '/data/EDB/MEG_AXCPT_Flanker/data/'; %top level directory for the BIDS output
    cfg.sub                     = sdan; %subject name
    cfg.run                     = ii; %optional
%     cfg.ses                     = '02'; %optional
    cfg.datatype                = 'meg'; %can be any of 'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'
    % cfg.participants.age        = meginfo.Age(iiS);
    % cfg.participants.sex        = meginfo.Sex(iiS); %'m' or 'f'
    % cfg.participants.group      = meginfo.Group(iiS);

    
    
    if strcmp(task,'flanker')
        cfg.method                  = 'copy';
        cfg.outputfile              = filenames{ii};
        
        cfg.task                    = 'flanker'; %task name is required for functional data
        cfg.TaskName                = cfg.task;
        cfg.TaskDescription         = 'Flanker task, 4 blocks of 48 trials each, congruent/incongruent, left and right stimuli with equal likelihood.';
    else
%         cfg.method                  = 'copy';
%         cfg.outputfile              = filenames{ii};
        cfg.method                  = 'decorate';
        cfg.outputfile              = filenames_old{ii};
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

%%
for iiS = 48
sdan = num2str(meginfo.SDAN(iiS));

cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan])
if ~exist('anat','dir')
    mkdir anat 
    mkdir meg
    !mv *.ds meg
    !mv *.nii *.json anat
end

cd anat
% convert anatomical to bids
mrNames = dir;
mrNames = {mrNames.name};
ind = contains(mrNames,'.nii');
if nnz(ind)>1
    keyboard %dbcont
end
% Easier to rename than using datas2bids
% cfg =[];
% cfg.dataset                 = mrNames{ind};
% cfg.sub                     = sdan;
% cfg.method                  = 'copy';
% cfg.bidsroot                = '/data/MBDU/bids/meg_mmi3/';
% cfg.datatype                = 'T1w';
% cfg.acq                     = 'mprage';
% cfg.writejson               = 'no';
% data2bids(cfg);

indEx = strfind(mrNames{ind},'.');
mrEx = mrNames{ind}(indEx(1):end);  % file extension (.nii or .nii.gz)
eval(sprintf('!mv %s sub-%s_acq-mprage_T1w%s',mrNames{ind},sdan,mrEx))

ind = contains(mrNames,'.json');
eval(sprintf('!mv %s sub-%s_acq-mprage_T1w.json',mrNames{ind},sdan))
end