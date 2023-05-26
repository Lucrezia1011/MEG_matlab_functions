addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

meginfo = readtable('~/meg-mmi/meg-mmi_main/MEG_participantsinfo.csv');

%% Move fiducial files to derivatives folder
for iiS = 1:size(meginfo,1)
    sdan = num2str(meginfo.SDAN(iiS));

    cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])
    mrNames = dir;
    mrNames = {mrNames.name};
    ind = find(contains(mrNames,'.tag'));
    if nnz(ind)==1
        eval(['!mv ', mrNames{ind} ,' /data/MBDU/MEG_MMI3/data/derivatives/sub-' sdan '/'])
    end
end

%%
delete('/data/MBDU/bids/meg_mmi3/dataset_description.json')
clc
for iiS = 1:size(meginfo,1)
sdan = num2str(meginfo.SDAN(iiS));

cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan])
delete('*.tsv' )

cd meg
% bids naming convention
% MEG data:  sub-<label>_task-<label>_run-<index>_meg.ds
% MRI data:  sub-<label>_T1w
% behavioral :  sub-<label>_task-<label>_run-<index>_events.tsv

mrNames = dir;
mrNames = {mrNames.name};
ind = find(contains(mrNames,'.ds'));
%%

for nInd = 1:length(ind)
cfg=[];

cfg.method                  = 'decorate';
cfg.meg.writesidecar        ='replace';
cfg.coordsystem.writesidecar = 'replace';
cfg.channels.writesidecar  ='replace';
cfg.events.writesidecar  = 'replace';
% cfg.writejson    = 'replace';%string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
% cfg.writetsv     = 'replace';%string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
% cfg.proc                    = 'raw';
cfg.dataset                 = mrNames{ind(nInd)};
cfg.bidsroot                = '/data/MBDU/bids/meg_mmi3/'; %top level directory for the BIDS output
cfg.sub                     = sdan; %subject name
cfg.run                     = str2double(cfg.dataset(strfind(cfg.dataset ,'run')+4)); %optional
cfg.datatype                = 'meg'; %can be any of 'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'
cfg.participants.age        = meginfo.Age(iiS);
cfg.participants.sex        = meginfo.Sex(iiS); %'m' or 'f'
cfg.participants.group      = meginfo.Group(iiS);

cfg.InstitutionName             = 'National Istitutes of Health';
cfg.InstitutionAddress          = '9000 Rockville Pike, Bethesda, MD, USA';
cfg.InstitutionalDepartmentName = 'NIH MEG core ';
% cfg.Manufacturer                = 'CTF Systems, Inc., Canada';
% cfg.ManufacturersModelName      = '';
% cfg.DeviceSerialNumber          = '';
% cfg.SoftwareVersions            = '';

cfg.meg.SamplingFrequency             = 1200; % REQUIRED. Sampling frequency (in Hz) of all the data in the recording, regardless of their type (e.g., 2400)
cfg.meg.PowerLineFrequency            = 60; % REQUIRED. Frequency (in Hz) of the power grid at the geographical location of the MEG instrument (i.e. 50 or 60)
cfg.meg.DewarPosition                 = 'upright'; % REQUIRED. Position of the dewar during the MEG scan: "upright", "supine" or "degrees" of angle from vertical: for example on CTF systems, upright=15??, supine = 90??.
cfg.meg.SoftwareFilters               = [];
cfg.meg.SoftwareFilters.SpatialCompensation.GradientOrder = "^3"; % REQUIRED. List of temporal and/or spatial software filters applied, orideally key:valuepairsofpre-appliedsoftwarefiltersandtheir parameter values: e.g., {"SSS": {"frame": "head", "badlimit": 7}}, {"SpatialCompensation": {"GradientOrder": Order of the gradient compensation}}. Write "n/a" if no software filters applied.
cfg.meg.DigitizedLandmarks            = true; % REQUIRED. Boolean ("true" or "false") value indicating whether anatomical landmark points (i.e. fiducials) are contained within this recording.
cfg.meg.DigitizedHeadPoints           = false; % REQUIRED. Boolean ("true" or "false") value indicating whether head points outlining the scalp/

cfg.dataset_description.writesidhecar        = 'yes';
cfg.dataset_description.Name                = 'MEG closed loop mood induction with monetary gambling in adolescents';
cfg.dataset_description.BIDSVersion         = '1.2.0'; %?
cfg.dataset_description.License             = 'CC0';
cfg.dataset_description.Authors             = {"Lucrezia Liuzzi", ...
    "Emilio Alejandro Valadez","Andre Zugman","Daniel Pine"}; %string or cell-array of strings
cfg.dataset_description.Acknowledgements    = ...
    ['This research was supported by the Intramural Research Program of the ',...
    'National Institute of Mental Health National Institutes of Health (NIH).',...
    ' This work used the computational resources of the NIH HPC ',...
    '(high-performance computing) Biowulf cluster (http://hpc.nih.gov).'];
cfg.dataset_description.HowToAcknowledge    = ['Please cite: '];
cfg.dataset_description.Funding             = {['Intramural Research Program of ',...
    'the National Institute of Mental Health National Institutes of Health (NIH) (Grant No. ZIA-?? [to DP])']};
cfg.dataset_description.ReferencesAndLinks  = {'https://www.biorxiv.org/content/10.1101/2021.03.04.433969v1'};
cfg.dataset_description.DatasetDOI          = 'doi: 10.1101/2021.03.04.433969';

cfg.outputfile              = cfg.dataset;
% if nInd == 1
%     cfg.task                    = 'mmi3'; %task name is required for functional data
%     cfg.TaskName                = cfg.task;
%     cfg.TaskDescription         = 'Mood Machine Interface, 3-levels: high-low-high mood targets';
% else
%     cfg.task                    = 'rest';%'mmi3'; %task name is required for functional data
%     cfg.TaskName                = cfg.task;
%     cfg.TaskDescription         = 'Eyes open resting state scan.';
% end


cfg.task                    = cfg.dataset(strfind(cfg.dataset ,'task')+5:strfind(cfg.dataset ,'run')-2);%task name is required for functional data
cfg.TaskName                = cfg.task;
if strcmp(cfg.task,'flanker')
    cfg.TaskDescription         = 'Flanker';
else
    cfg.TaskDescription         = 'AX-CPT';
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
% eval(sprintf(['!newDs -includeBad -includeBadChannels -includeBadSegments ',...
%     '%s /data/MBDU/bids/meg_mmi3/sub-%s/meg/%s'],mrNames{ind(1)},sdan))

%%

sdan = num2str(meginfo.SDAN(iiS));

cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/anat'])

% convert anatomical to bids
mrNames = dir;
mrNames = {mrNames.name};
ind = contains(mrNames,'.nii');


if nnz(ind)>0 % check if there is an anatomical file
 
    if ~contains(mrNames{ind},'.gz') && nnz(ind)>0
        eval(['!gzip ',mrNames{ind}])
        mrNames = dir;
        mrNames = {mrNames.name};
        ind = contains(mrNames,'.nii');
    end
    
    cfg=[];
    
    indj = contains(mrNames,'.json');
    json_name = mrNames{indj};
    ft_info('reading %s\n', json_name);
    ft_hastoolbox('jsonlab', 1);
    json = loadjson(json_name);
    json = ft_struct2char(json); %
    cfg.mri = json;
    
    fids_name = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sdan,'/sub-',sdan,'_fiducials.tag'];

    if exist(fids_name,'file')
        fileID = fopen(fids_name,'r');
        fids_char = fscanf(fileID,'%c');
        fclose(fileID);

        fids_inds = zeros(3,3);

        % 4 lines of 66 characters each
        for iF = 1:3
            fids_inds(iF,1) = str2double(fids_char(66*iF+(18:28))); 
            fids_inds(iF,2) = str2double(fids_char(66*iF+(30:40)));
            fids_inds(iF,3) = str2double(fids_char(66*iF+(42:52)));
        end
        cfg.mri.AnatomicalLandmarkCoordinates = {'NAS:',fids_inds(1,:),'LPA:',fids_inds(2,:),'RPA:',fids_inds(3,:)};
        cfg.mri.AnatomicalLandmarkCoordinateUnits = 'mm';
        cfg.mri.AnatomicalLandmarkCoordinateSystem = 'LPS';
        cfg.mri.AnatomicalLandmarkCoordinateSystemDescription = 'Anatomical Landmarks registered with AFNI software Edit Tagset plugin. Defaults to LPS coordinate system.';

    end

    cfg.method                  = 'decorate';
    cfg.mri.writesidecar        = 'replace';

    cfg.writejson    = 'replace';%string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
    % cfg.writetsv     = 'replace';%string, 'yes', 'replace', 'merge' or 'no' (default = 'yes')
    % cfg.proc                    = 'raw';
    cfg.dataset                 = mrNames{ind};
    cfg.bidsroot                = '/data/MBDU/MEG_MMI3/data/bids/'; %top level directory for the BIDS output
    cfg.sub                     = sdan; %subject name
    cfg.run                     = 1; %optional
    cfg.datatype                = 'anat'; %can be any of 'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'
    cfg.participants.age        = meginfo.Age(iiS);
    cfg.participants.sex        = meginfo.Sex(iiS); %'m' or 'f'
    cfg.participants.group      = meginfo.Group(iiS);
    cfg.outputfile              = cfg.dataset;

   
    cfg.scans.acq_time          = json.AcquisitionDateTime;
    
    data2bids(cfg);
    % Make correct symbolic link to mri json file
%     delete( ['/data/MBDU/MEG_MMI3/data/bids_defaced/sub-' sdan '/anat/' json_name])
%     eval(['!ln -s /data/MBDU/bids/meg_mmi3/sub-' sdan '/anat/' json_name ...
%         ' /data/MBDU/MEG_MMI3/data/bids_defaced/sub-' sdan '/anat/' json_name  ])
end
end



%% Find empty files in .ds directories
filedels2 = [];
n = 0;
for iiS = 1:size(meginfo,1)
    sdan = num2str(meginfo.SDAN(iiS));

    dirname = ['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg/'];
    cd(dirname)
    mrNames = dir;
    mrNames = {mrNames.name};
    ind = find(contains(mrNames,'.ds'));
    %%

    for nInd = 1:length(ind)
        cd(mrNames{ind(nInd)})
        d = dir;
        d(1:2) = [];
        % find empty files
        b = [d.bytes];
        emptyFiles = {d(b == 0).name};
%         ii = find(contains(emptyFiles,'.bak')); % fidels 57 files: 54 files gave problems on openeuro
%         ii = find(contains(emptyFiles,'bad.segments')); fidels2 %
%         deleting bad.segments does not help
%         ii = find(~isfolder(emptyFiles) & ~contains(emptyFiles,'.eeg')) ;
        ii = find(~isfolder(emptyFiles) & contains(emptyFiles,'.eeg')) ; % 30 files
%         %fidels3 deleting badChannels does not help
        for iix = 1:length(ii)
            n = n+1;
            filedels2{n} = [dirname,mrNames{ind(nInd)},'/', emptyFiles{ii(iix)}];
            delete(filedels2{n})
        end
        cd ..
    end
end

save('~/openneuro_dels','filedels2','-append')

% save('~/openneuro_dels','filedels','-append')

%% Put back deleted empty files

load ~/openneuro_dels.mat

backuppath = '/data/MBDU/.snapshot/weekly.2021-03-02_0717/';

for n = 1:length(filedels)
    filebck = [backuppath, filedels{n}(12:end)];
    cm = ['!cp ',filebck,' ', filedels{n}];
    eval(cm)
end

for n = 1:length(filedels2)
    filebck = [backuppath, filedels2{n}(12:end)];
    cm = ['!cp ',filebck,' ', filedels2{n}];
    eval(cm)
end
% 
% for n = 1:length(filedels3)
%     filebck = [backuppath, filedels3{n}(12:end)];
%     cm = ['!cp ',filebck,' ', filedels3{n}];
%     eval(cm)
% end