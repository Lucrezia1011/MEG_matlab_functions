%function AXCPT_process_source(filenames,datapath,processingfolder,bvfiles,fids_name)

addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

sub = '24590';

datapath = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/meg/'];

processingfolder = ['/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
if ~exist(processingfolder,'dir')
    mkdir(processingfolder)
end

filenames = cell(1,2);
for ii = 1:2
    filenames{ii} = ['sub-',sub,'_task-axcpt_run-',num2str(ii),'_meg.ds'];
end
mri_name = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/anat/sub-',sub,'_acq-mprage_T1w.nii'];
fids_name = ['sub-',sub,'_fiducials.tag'];

close all

%% Empty rooms


cd(datapath)
cd(filenames{1})
fileID = fopen([filenames{1}(1:end-2),'infods'],'r');
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

d = dir('/data/liuzzil2/MEG_AXCPT_Flanker/data/emptyroom/');
emptyroom = []; jj = 0;
for ii = 3:length(d)
    if contains(d(ii).name, TaskDate(1:8))
        jj = jj + 1;
        emptyroom{jj} = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/emptyroom/',d(ii).name];
    end
end

 
highpass = 0.5;
lowpass = 120;
icaopt = 1;
plotopt = 0;

noiseC = zeros(271,271,length(emptyroom));
for ii = 1:length(emptyroom)
    cfg = [];
    cfg.dataset = emptyroom{ii};
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';
    cfg.demean = 'yes';
    cfg.detrend = 'no';
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [lowpass, highpass]; % With notch filter 60Hz
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [58 62; 118 122; 178 182]; % With notch filter 60Hz

    data_empty = ft_preprocessing(cfg);
    emptyC = zeros(271,271,length(data_empty.trial));
    for t = 1:length(data_empty.trial)
        emptyC = cov(data_empty.trial{t}');
    end
    noiseC(:,:,ii) = mean(emptyC,3);
end
noiseC = mean(noiseC,3);


%% Standard pre-processing

cd(datapath)


resultsfolder = [processingfolder,'axcpt/'];

for ii = 1:length(filenames)
    filename = filenames{ii};
    sub = filename(5:9);
    
    [data,BadSamples] = preproc_bids(filename,processingfolder,highpass,lowpass,icaopt,plotopt);
    % eyelink
    if any(strcmp(data.hdr.label,'UADC009'))
        cfg = [];
        cfg.dataset = filename;
        cfg.continuous = 'yes';
        cfg.channel = {'UADC009';'UADC010';'UADC013'};
        eyed = ft_preprocessing(cfg);
        eyelink = eyed.trial{1};
        eyelink(:,BadSamples) = [];
    end
     
    % Read events
    [samples_T,trig_sample,buttonpress] = matchTriggers( filename, BadSamples);
    
    trig_sample.sample = samples_T(trig_sample.sample);
    trig_sample.type(trig_sample.sample == 0) = [];
    trig_sample.value(trig_sample.sample == 0) = [];
    trig_sample.sample(trig_sample.sample == 0) = [];
    
    buttonpress.left = samples_T(buttonpress.UADC006);
    buttonpress.left(buttonpress.left == 0) = [];
    buttonpress.right = samples_T(buttonpress.UADC007);
    buttonpress.right(buttonpress.right == 0) = [];
    
    if data.hdr.nSamplesPre == 0
        time=  1:length(samples_T);
        time(samples_T == 0) =[];
    else
        time= repmat( ((1:data.hdr.nSamples) - data.hdr.nSamplesPre),[1,data.hdr.nTrials] );
        
    end
  
    
    if ii == 1
        trigall = trig_sample;
        buttonresp = buttonpress;
        dataall = data;
        timeall = time;
        if any(strcmp(data.hdr.label,'UADC009'))
            eyeall = eyelink;
        end
          
        mri_name = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-01/anat/sub-',sub,'_acq-mprage_T1w.nii'];
        if ~exist(mri_name,'file')
            mri_name = [mri_name,'.gz'];
        end
        fids_file =  [datapath(1:end-4),'anat/',fids_name];
        mri = fids2ctf(mri_name,fids_file,0);
        
        mneres =5124; % 5124, 8196, 20484
        gridl =mniLeadfields_singleshell(filename,processingfolder,mneres,mri); % calculate leadfields on MNI grid
        load([processingfolder,'/headmodel_singleshell.mat']);
        
    else
        
        if any(strcmp(data.hdr.label,'UADC009'))
            eyeall = cat(2,eyeall,eyelink);
        end
        
        % realign sensors to first recording
        cfg = [];
        cfg.template          = dataall.grad;
        cfg.headmodel  =  vol;
        cfg.inwardshift = 1; 
        [data] = ft_megrealign(cfg, data);
  
        
        trigall.sample = cat(1, trigall.sample, trig_sample.sample + dataall.sampleinfo(2));
        trigall.value = cat(1, trigall.value, trig_sample.value);
        trigall.type = cat(1, trigall.type, trig_sample.type);
        
        buttonresp.right = cat(1, buttonresp.right, buttonpress.right + dataall.sampleinfo(2));
        buttonresp.left = cat(1, buttonresp.left, buttonpress.left + dataall.sampleinfo(2));
        
        dataall.trial{1} = cat(2,dataall.trial{1},data.trial{1}(1:length(dataall.label),:));
        dataall.time{1} = cat(2,dataall.time{1},data.time{1} + dataall.time{1}(end) + 5);
        dataall.sampleinfo(2) = length(dataall.time{1});
        timeall = cat(2, timeall, time);
        
    end
    
end

data = dataall;
close all
%% Find samples
cue = 'A';
% buttonpress to cue
cuesampA = trigall.sample(strcmp(trigall.type, cue) ) ;

[buttonsampA,rtA] = match_responses(cuesampA, buttonresp, 'left', data.fsample);

cue = 'B';
% buttonpress to cue
cuesampB = trigall.sample(strcmp(trigall.type, cue) ) ;
[buttonsampB,rtB] = match_responses(cuesampB, buttonresp, 'left', data.fsample);

probe = 'AY';
%buttonpress to probe
probesampAY = trigall.sample(strcmp(trigall.type, probe)) ;

buttonsampAYcomm = match_responses(probesampAY, buttonresp, 'right', data.fsample);
buttonsampAYcorr = match_responses(probesampAY, buttonresp, 'left', data.fsample);

probe = 'AX';
%buttonpress to probe
probesampAX = trigall.sample(strcmp(trigall.type, probe)) ;

buttonsampAXcorr = match_responses(probesampAX, buttonresp, 'right', data.fsample);
buttonsampAXcomm = match_responses(probesampAX, buttonresp, 'left', data.fsample);


probe = 'BX';
%buttonpress to probe
probesampBX = trigall.sample(strcmp(trigall.type, probe)) ;

buttonsampBXcomm = match_responses(probesampBX, buttonresp, 'right', data.fsample);
buttonsampBXcorr = match_responses(probesampBX, buttonresp, 'left', data.fsample);
%% Filters
freq = [120];
filt_order = []; % default

if size(freq,2) == 1
    data_filt = ft_preproc_lowpassfilter(dataall.trial{1}, data.fsample,freq,filt_order,'but');
    data.trial{1} = data_filt;
    clear data_filt
elseif size(freq,2) == 2
    data_filt = ft_preproc_bandpassfilter(dataall.trial{1}, dataall.fsample,freq,filt_order,'but');
    data.trial{1} = data_filt;
    clear data_filt
end
 
    %% Define trials
        
    twind = [-0.5 1];
   
 
    [datacueAresp,ttdel] = define_trials(buttonsampA,data,timeall,twind,1);

  
    [dataprobeAY,ttdel] = define_trials(probesampAY ,data,timeall,twind,1);
    

    [datacueA,ttdel] = define_trials(cuesampA,data,timeall,twind,1);
    [datacueB,ttdel] = define_trials(cuesampB,data,timeall,twind,1);
    [datacue,ttdel] = define_trials([cuesampA;cuesampB],data,timeall,twind,1);
    
    [dataprobeAX,ttdel] = define_trials(probesampAX,data,timeall,twind,1);
    [dataprobeBX,ttdel] = define_trials(probesampBX,data,timeall,twind,1);
    
    [dataproberespComm,ttdel] = define_trials([buttonsampAYcomm;buttonsampAXcomm],data,timeall,twind,1);
    [dataproberespAYComm,ttdel] = define_trials([buttonsampAYcomm],data,timeall,twind,1);
    
    [dataproberespAXCorr,ttdel] = define_trials(buttonsampAXcorr,data,timeall,twind,1);
    [dataproberespAYCorr,ttdel] = define_trials(buttonsampAYcorr,data,timeall,twind,1);
    [dataprobeArespCorr,~] =  define_trials([buttonsampAXcorr;buttonsampAYcorr],data,timeall,twind,1);
    
    cfg = [];
    cfg.demean        = 'yes';
    cfg.baselinewindow = [-inf, -0];
    
    if strcmp(cfg.demean, 'yes')
        datacueA = ft_preprocessing(cfg,datacueA);
        datacueB = ft_preprocessing(cfg,datacueB);
        datacue = ft_preprocessing(cfg,datacue);
        dataprobeAX = ft_preprocessing(cfg,dataprobeAX);
        dataprobeAY = ft_preprocessing(cfg,dataprobeAY);
        dataprobeBX = ft_preprocessing(cfg,dataprobeBX);
        
        dataproberespComm = ft_preprocessing(cfg,dataproberespComm);
        dataproberespAYComm = ft_preprocessing(cfg,dataproberespAYComm);
        
        dataproberespAXCorr = ft_preprocessing(cfg,dataproberespAXCorr);
        dataproberespAYCorr = ft_preprocessing(cfg,dataproberespAYCorr);
    end
    

   
    %% Co-register MRI
    
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = [-inf 0]; %it will calculate the covariance matrix
    cfg.keeptrials    = 'no';
    % on the timepoints that are
    % before the zero-time point in the trials
    cfg.trials = 1:length(datacueB.trial);
    tlckA = ft_timelockanalysis(cfg, datacueA);
    tlckB = ft_timelockanalysis(cfg, datacueB);
    
    cfg.trials = 1:length(dataprobeAY.trial);
    tlckAX = ft_timelockanalysis(cfg, dataprobeAX);
    tlckAY = ft_timelockanalysis(cfg, dataprobeAY);
    
    cfg               = [];
    cfg.method        = 'mne';
    cfg.sourcemodel   = gridl;
    cfg.headmodel     = vol;
    cfg.mne.prewhiten = 'yes';
    cfg.mne.lambda    = 3;
    cfg.mne.scalesourcecov = 'yes';
    cfg.mne.keepfilter = 'no';

    sourceA          = ft_sourceanalysis(cfg,tlckA);
    sourceB         = ft_sourceanalysis(cfg, tlckB);
    sourceAX          = ft_sourceanalysis(cfg,tlckAX);
    sourceAY          = ft_sourceanalysis(cfg,tlckAY);

%% plot MNE

cfg = [];
cfg.projectmom = 'yes';
sdA  = ft_sourcedescriptives(cfg,sourceA);
sdB = ft_sourcedescriptives(cfg, sourceB);

sdAX  = ft_sourcedescriptives(cfg,sourceAX);
sdAY  = ft_sourcedescriptives(cfg,sourceAY);


% 500 ms after the zero time-point
sdDIFF         = sdAX;
sdDIFF.avg.pow = sdAY.avg.pow - sdAX.avg.pow;

for n = 1:size(sdDIFF.avg.pow,1)
   sdDIFF.avg.pow(n,:) = smooth(sdDIFF.avg.pow(n,:),10) ;
end

powdist = sdDIFF.avg.pow(:,sdDIFF.time < 0);
min(powdist(:))*2
max(powdist(:))*2


cfg = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,sdDIFF);


d = sdDIFF.avg.pow;
for n = 1:size(d,1)
   d(n,:) = smooth(d(n,:),10) ;
end
[coeff,score]= pca(d');

sourceAr = sourceA;
sourceAr.avg.pow = coeff;
for jj = 1:10
m=sourceAr.avg.pow(:, jj); % plotting the result at the 450th time-point that is
figure; set(gcf,'position',[182   539   888   363],'color','w')
subplot(1,2,1)
ft_plot_mesh(sourceAr, 'vertexcolor', m); %caxis([0 1]*4e-22); 
caxis([-1 1]*0.05); 
view([150 20]); h = light; set(h, 'position', [0 1 0.2]); lighting gouraud; material dull
subplot(1,2,2)
plot(sourceAr.time, score(:,jj)); grid on
end


%% Beamfomer leadfields
% Co-register MRI

    gridres = 5; % 5mm grid
    gridl =mniLeadfields_multiSpheres(filenames{ii},processingfolder,gridres,mri); % calculate leadfields on MNI grid

    % For power
        % Load standard brain for plotting
    %     mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
    mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates
    mri_mni.coordsys = 'mni';
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    sourcemodel.coordsys = 'mni';

%% Beamfomer
    if isfield(data.cfg,'component')
        icacomps = length(data.cfg.component);
    else
        icacomps = 0;
    end
    C = cov(data.trial{1}');

    nchans = length(data.label);
    % E = svd(C);
    % noiseSVD = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
    E = svd(noiseC);
    noiseSVD = eye(nchans)*E(end);
    mu  =4;
    Cr = C + mu*noiseSVD; % old normalization
    %     Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value

    L = gridl.leadfield(gridl.inside);

    weigths_file = sprintf('%s%s/weights_multiSpheres_%dmm_regmu%d.mat',processingfolder,filename(1:end-3),gridres,mu);
    if ~exist(weigths_file,'file')
        W = cell(size(L));
        Wdc = cell(size(L));
        parfor l = 1:length(L)

            lf = L{l}; % Unit 1Am

            % G O'Neill method, equivalent to fieldtrip
            [v,d] = svd(lf'/Cr*lf);
            d = diag(d);
            jj = 2;

            lfo = lf*v(:,jj); % Lead field with selected orientation

            w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
            Wdc{l} = w;
        % no depth correction as we later divide by noise
            w = Cr\lfo / (lfo'/Cr*lfo) ; % weights
            W{l} = w;

            if mod(l,300) == 0
                clc
                fprintf('SAM running %.1f\n',...
                    l/length(L)*100)
            end

        end
        save(weigths_file,'W','Wdc')
    else
        load(weigths_file)
    end



%% Beamfomer weights

% tlckA = ft_timelockanalysis([], dataprobeAX);
% tlckB = ft_timelockanalysis([], dataprobeAY);
% figure; subplot(311); plot(tlckA.time, abs(tlckA.avg) )
%  subplot(312); plot(tlckA.time, abs(tlckB.avg) )
%  subplot(313); plot(tlckA.time, abs(tlckA.avg - tlckB.avg) )

% tt= 0.41; 
% [~,iit] = min(abs(tlckA.time -tt));
% probeAX =  tlckA.avg(:,iit );
% probeAY=  tlckB.avg(:,iit);

timew = [0.1 0.6];
dataA = define_trials(probesampAY,data,timeall,timew,1);
dataC = define_trials(probesampAX,data,timeall,timew,1);

% dataA = define_trials(cuesampB,data,timeall,timew,1);
% dataC = define_trials(cuesampA,data,timeall,timew,1);
% 
% 
Ca = cov(cell2mat(dataA.trial)');
% Cc = cov(cell2mat(dataC.trial)');
Cc = cov(cell2mat(dataC.trial(randperm(length(dataC.trial),length(dataA.trial))))');


P = cell(size(L));
Pc = cell(size(L));


for ii = 1:length(L)
    
    w = W{ii};
    
%         P{ii} = (w'*probeAY); % needs depth corrected weights
%         Pc{ii} = (w'*probeAX);
   
    
        P{ii} = ((w'*Ca*w) - (w'*Cc*w) ) / (2*w'*noiseC*w);
%         P{ii} = (w'*Ca*w)  / (w'*noiseC*w);
%         Pc{ii} = (w'*Cc*w)  / (w'*noiseC*w);
    
    
end

P  = cell2mat(P)';
Pc  = cell2mat(Pc)';
%%
mriplot = 'mni';
plotopt = 'ortho';  % 'slice'; %'ortho'

crang = [-5 5];


T = zeros(gridl.dim);
T(gridl.inside) =(P);
sourceant =[];
sourceant.pow = T;
if strcmp(mriplot,'mni')
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    sourceout_Int.coordsys = 'mni';
    
else
    sourceant.dim = gridl.dim;
    sourceant.inside = gridl.inside;
    sourceant.pos = gridl.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    
end

% crang = max(sourceout_Int.pow(:))* [.3 1]; % 50 % of maximum

cfg = [];
cfg.method        = plotopt;
cfg.location   = [-12 23 56];

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
if strcmp(mriplot,'mni')
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
end
ft_sourceplot(cfg, sourceout_Int);
set(gcf,'position',[704   117   894   675],'color','w')


% title(sprintf('Probe AY - AX, %dms',round(tt*1e3)))

% set(gcf,'name',sprintf('CueA ERF peak at %dms ',round(tt*1e3)))
% saveas(gcf,sprintf('%s/ProbeAY-AX_ERF_%dms_ortho.jpg',resultsfolder,round(tt*1e3)))

% 
saveas(gcf,sprintf('%s/CueB-CueA_%d-%dms_%d-%dHz.jpg',resultsfolder,...
    timew(1)*1000, timew(2)*1000,freq(1),freq(2)))

% saveas(gcf,sprintf('%s/ProbeAY-AX_%d-%dms_%d-%dHz.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000,freq(1),freq(2)))

%% Beamformer fieldtrip

addpath /home/liuzzil2/fieldtrip-20190812/fieldtrip_private/

gridres = 5;
gridl =mniLeadfields_multiSpheres(filenames{1},processingfolder,gridres,mri); % calculate leadfields on MNI grid
load([processingfolder,'/headmodel_multiSpheres.mat']);

twind = [-0.5 1];
[datacue,~] = define_trials([cuesampA;cuesampB],data,timeall,twind,1);

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = twind; %it will calculate the covariance matrix
cfg.keeptrials    = 'yes';
% on the timepoints that are
% before the zero-time point in the trials
tlck = ft_timelockanalysis(cfg, datacue);
cfg = [];
cfg.method = 'lcmv';
cfg.grid = gridl;
cfg.headmodel = vol;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.lambda = '5%';
cfg.lcmv.kappa = length(data.label) - 4; % deleted ICA components
cfg.channel = {'MEG'};
cfg.senstype = 'MEG';
sourceavg = ft_sourceanalysis(cfg, tlck);

twind = [0.4 0.7];
[datawind,~] = define_trials([cuesampA;cuesampB],data,timeall,twind,1);
cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes'; % This is important as the single-trial source activity is calculated based on the covariance matrix
avgwind = ft_timelockanalysis(cfg,datawind);

% call to ft_sourceanalysis now applying the precomputed filters to pre and %post intervals
cfg = [];
cfg.method = 'lcmv';
cfg.sourcemodel = gridl;
cfg.sourcemodel.filter = sourceavg.avg.filter;
cfg.headmodel = vol;
cfg.keeptrials = 'yes';
cfg.rawtrial=  'yes';
% cfg.rawtrial      =  'yes';
sourcewind = ft_sourceanalysis(cfg, avgwind);

% Only works for volume beamformer
design = [ones(size(cuesampA)); ones(size(cuesampB))*2]';
cfg = [];
% cfg.dim         = source.dim;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_indepsamplesT';
cfg.parameter   = 'pow';
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000;
cfg.alpha       = 0.05; % note that this only implies single-sided testing
cfg.tail        = 0;
cfg.design(1,:) = [1:length(find(design==1)) 1:length(find(design==2))];
cfg.design(2,:) = design;
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: trials)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfg, sourcewind);


%% Beamformer fieldtrip plot

twind = [0 1];
[datapre,ttdel] = define_trials(cuesampB,data,timeall,[-diff(twind), 0],1);
[datapost,ttdel] = define_trials(cuesampB,data,timeall,twind,1);
cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes'; % This is important as the single-trial source activity is calculated based on the covariance matrix
avgpre = ft_timelockanalysis(cfg,datapre);
avgpst = ft_timelockanalysis(cfg,datapost);

% call to ft_sourceanalysis now applying the precomputed filters to pre and %post intervals
cfg = [];
cfg.method = 'lcmv';
cfg.grid = gridl;
cfg.sourcemodel.filter = sourceavg.avg.filter;
cfg.headmodel = vol;
sourcepre = ft_sourceanalysis(cfg, avgpre);
sourcepst = ft_sourceanalysis(cfg, avgpst);



% Plot Beamformer

T=sourcepst;
T.avg.pow=(sourcepst.avg.pow-sourcepre.avg.pow)./sourcepre.avg.pow;

cfg = [];
cfg.parameter = 'pow';
source_int  = ft_sourceinterpolate(cfg, T , mri);
source_int.pow(~source_int.inside) = 0;
% sourceout_Int.coordsys = 'mni';

source_int.mask = source_int.pow > max(source_int.pow(:))*.3; % 50 % of maximum
% crang = [];
crang = [-1 1]*max(source_int.pow(:));
cfg = [];
cfg.method        = 'ortho'; %'ortho'
cfg.location   = 'max';

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, source_int);
set(gcf,'position',[704   117   894   675],'color','w')
title(sprintf('ERF power at %d-%dms',twind(1)*1e3,twind(2)*1e3))
% title(sprintf('theta power %d-%dms after AX probe',timew(1)*1000, timew(2)*1000))
% saveas(gcf,sprintf('/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-24531/derivatives/theta_power_%d-%dms_AXprobe.jpg',...
%     timew(1)*1000, timew(2)*1000))

