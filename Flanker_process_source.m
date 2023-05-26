function Flanker_process(filenames,datapath,processingfolder,bvfiles,fids_name)

addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

sub = '24563';

datapath = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/meg/'];

processingfolder = ['/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
if ~exist(processingfolder,'dir')
    mkdir(processingfolder)
end
% fids_name = 'sub-24531_fiducials.tag';
% 

filenames = cell(1,2);
for ii = 1:2
    filenames{ii} = ['sub-',sub,'_task-flanker_run-',num2str(ii),'_meg.ds'];
end

fids_name = ['sub-',sub,'_fiducials.tag'];
mri_name = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/anat/sub-',sub,'_acq-mprage_T1w.nii'];

%% Standard pre-processing

cd(datapath)
 
highpass = 0.5;
lowpass = 120;
icaopt = 1;
plotopt = 0;

resultsfolder = [processingfolder,'flanker/'];


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

    % Read behavioral file
%     M = readtable(['../behavioral/',bvfiles{ii}]) ;
  
    % Read events
    [samples_T,trig_sample,buttonpress] = matchTriggers( filename, BadSamples);
    
    trig_sample.sample = samples_T(trig_sample.sample);
    trig_sample.sample(trig_sample.sample == 0) = [];
    if ~isempty(buttonpress)
        buttonpress.left = samples_T(buttonpress.UADC006);
        buttonpress.left(buttonpress.left == 0) = [];
        buttonpress.right = samples_T(buttonpress.UADC007);
        buttonpress.right(buttonpress.right == 0) = [];
    else
        ttype = M.Trialtype;
        ttypemeg = trig_sample.type(~strcmp(trig_sample.type, 'cue'));   
        zz = 1;
        mkeep = false(size(ttype));
        for jj = 1:length(ttype)-1
           mkeep(jj) = strcmpi(ttype{jj}, ttypemeg{zz}(1:2));
           zz = zz + mkeep(jj);
        end
        M= M(mkeep,:);
        flsamp = trig_sample.sample(~strcmp(trig_sample.type,'cue'));
        for cc = 1:size(M,1)
            if ~isempty(M.key_resp_rt{cc})
                M.key_resp_rt{cc} = round(str2double( M.key_resp_rt{cc}(2:end-1))*data.fsample) + flsamp(cc);
            end
        end
        
        buttonpress.left =( M.key_resp_rt(~cellfun(@isempty, (regexp(M.key_resp_keys,'3')))) );       
        buttonpress.left = cell2mat(buttonpress.left);
        
        buttonpress.right =( M.key_resp_rt(~cellfun(@isempty, (regexp(M.key_resp_keys,'1')))) );
        buttonpress.right = cell2mat(buttonpress.right);
        
        
    end
    if data.hdr.nSamplesPre == 0
        time= samples_T';
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
        
%         mri_name = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-01/anat/sub-',sub,'_acq-mprage_T1w.nii'];
       
        if ~exist(mri_name,'file')
            mri_name = [mri_name,'.gz'];
        end
        fids_file =  [datapath(1:end-4),'anat/',fids_name];
        mri = fids2ctf(mri_name,fids_file,0);
        
        mneres =20484; % 5124, 8196, 20484
        gridl =mniLeadfields_singleshell(filename,processingfolder,mneres,mri); % calculate leadfields on MNI grid
        load([processingfolder,'/headmodel_singleshell.mat']);
        
        
    else
        
        if any(strcmp(data.hdr.label,'UADC009'))
            eyeall = cat(2,eyeall,eyelink);
        end
        % realign sensors to first recording
%         cfg = [];
%         cfg.template          = dataall.grad;
%         cfg.headmodel  =  vol;
%         cfg.inwardshift = 1; 
%         [data] = ft_megrealign(cfg, data);
          
        
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


%% 
filt_order = []; % default
data_filt = ft_preproc_lowpassfilter(dataall.trial{1}, dataall.fsample,50,filt_order,'but');
% data_filt = ft_preproc_bandpassfilter(dataall.trial{1}, dataall.fsample,[13 30],filt_order,'but');
data = dataall;
data.trial{1} = data_filt;
clear data_filt

%% TFS

data = dataall;

timewresp = [-1 2];   

stimname = 'resp';
if strcmp(stimname,'resp' )
    dataIcorr = define_trials([respRIcorr;respLIcorr],data,timeall,timewresp,1);
    dataIcomm = define_trials([respRIcomm;respLIcomm],data,timeall,timewresp,1);
    dataCcorr = define_trials([respRCcorr;respLCcorr],data,timeall,timewresp,1);
else
    dataIcorr = define_trials([flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
    dataIcomm = define_trials([flanRIcomm;flanLIcomm],data,timeall,timewresp,1);
    dataCcorr = define_trials([flanRCcorr;flanLCcorr],data,timeall,timewresp,1);
    datacorr = define_trials([flanRCcorr;flanLCcorr;flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
end



cfg                 = [];
cfg.method          = 'mtmconvol';
cfg.taper           = 'hanning';
cfg.channel         = 'MEG';

% set the frequencies of interest
cfg.foi             = 1:2:100;
% set the timepoints of interest: from -0.8 to 1.1 in steps of 100ms
cfg.toi             = -1:0.1:2;

% set the time window for TFR analysis: constant length of 200ms
cfg.t_ftimwin       = 0.2 * ones(length(cfg.foi), 1);

% average over trials
cfg.keeptrials      = 'no';

% pad trials to integer number of seconds, this speeds up the analysis
% and results in a neatly spaced frequency axis
cfg.pad             = 4;
freqIcorr           = ft_freqanalysis(cfg, dataIcorr);
freqIcomm           = ft_freqanalysis(cfg, dataIcomm);
freqdiff = freqIcorr;
freqdiff.powspctrm = freqIcomm.powspctrm - freqIcorr.powspctrm;

figure; set(gcf,'color','w')
changroups = {'MLF';'MRF';'MLC';'MRC';'MLT';'MRT';'MLP';'MRP';'MLO';'MRO'};
cfg                 = [];
cfg.baseline        = [-1 -0.3];
cfg.baselinetype    = 'relchange';
cfg.zlim            = [-1 1];
cfg.ylim            = [0 50];
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,freqIcorr)
    title(changroups{n});
end

figure; set(gcf,'color','w')
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,freqIcomm)
    title(changroups{n});
end

figure; set(gcf,'color','w')
cfg                 = [];
cfg.zlim            ='maxabs';
cfg.ylim            = [0 50];
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,freqdiff)
    title(changroups{n});
end

% saveas(gcf,sprintf('%s/Corr_TFS_sens_%s.jpg',resultsfolder,stimname))
%%


timewresp = [-0.5 1];   
flan = 'RI';
resp = 'left';
% buttonpress to cue
cuesampr = trigall.sample(strcmp(trigall.type, flan) ) ;
[respRIcorr,flanRIcorr] = resp2trial(resp,cuesampr,buttonresp,data.fsample);

resp = 'right';
[respRIcomm,flanRIcomm] = resp2trial(resp,cuesampr,buttonresp,data.fsample);


flan = 'LI';
resp = 'right';
% buttonpress to cue
cuesampl = trigall.sample(strcmp(trigall.type, flan) ) ;
[respLIcorr,flanLIcorr] = resp2trial(resp,cuesampl,buttonresp,data.fsample);

resp = 'left';
[respLIcomm,flanLIcomm] = resp2trial(resp,cuesampl,buttonresp,data.fsample);


flan = 'LC';
resp = 'left';
% buttonpress to cue
cuesamp = trigall.sample(strcmp(trigall.type, flan) ) ;
[respLCcorr,flanLCcorr] = resp2trial(resp,cuesamp,buttonresp,data.fsample);

flan = 'RC';
resp = 'right';
% buttonpress to cue
cuesamp = trigall.sample(strcmp(trigall.type, flan) ) ;
[respRCcorr,flanRCcorr] = resp2trial(resp,cuesamp,buttonresp,data.fsample);

stimname = 'resp';
if strcmp(stimname,'resp' )
    dataIcorr = define_trials([respRIcorr;respLIcorr],data,timeall,timewresp,1);
    dataIcomm = define_trials([respRIcomm;respLIcomm],data,timeall,timewresp,1);
    dataCcorr = define_trials([respRCcorr;respLCcorr],data,timeall,timewresp,1);
else
    dataIcorr = define_trials([flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
    dataIcomm = define_trials([flanRIcomm;flanLIcomm],data,timeall,timewresp,1);
    dataCcorr = define_trials([flanRCcorr;flanLCcorr],data,timeall,timewresp,1);
    datacorr = define_trials([flanRCcorr;flanLCcorr;flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
end


% Ca = cov(cell2mat(dataI.trial)');
% Cc = cov(cell2mat(dataC.trial)');


%% Plot sens



% cfg = [];
% cfg.demean = 'yes';
% cfg.baselinewindow = [-0.1, -0.0];
% dataIcomm = ft_preprocessing(cfg,dataIcomm);
% dataIcorr = ft_preprocessing(cfg,dataIcorr);

cfg= [];
cfg.keeptrials         = 'no';
erpComm = ft_timelockanalysis(cfg, dataIcomm);
% erpComm.avg = squeeze(mean(erpComm.trial,1));
erpICorr = ft_timelockanalysis(cfg, dataIcorr);
% erpICorr.avg = squeeze(mean(erpICorr.trial,1));
erpCCorr = ft_timelockanalysis(cfg, dataCcorr);
% erpCCorr.avg = squeeze(mean(erpCCorr.trial,1));

if strcmp(stimname,'resp' )
    erpDiff = erpComm;
    erpDiff.avg = erpComm.avg - erpICorr.avg;
else
    erpDiff = erpICorr;
    erpDiff.avg = erpICorr.avg - erpCCorr.avg;
end
figure; set(gcf,'position',[644   457   998   433])
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.parameter = 'avg';
cfg.interpolatenan = 'no';
cfg.xlim  = [-0.05 0.];
cfg.zlim  = [-1 1]*1e-13;
cfg.colorbar  = 'yes';
subplot(122)
ft_topoplotER(cfg, erpComm)
subplot(121)
ft_topoplotER(cfg, erpICorr)

cfg = [];
cfg.latency = [0.100 0.110];
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.resolution = 1;
cfg.unit = 'cm';
cfg.gridsearch = 'yes';
cfg.headmodel = headmodel_meg;
cfg.senstype = 'meg';
cfg.channel = {'MEG*2', 'MEG*3'};
source_planar = ft_dipolefitting(cfg, timelock_all);
% 
figure
clf
changroups = {'MLF';'MRF';'MLC';'MRC';'MLT';'MRT';'MLP';'MRP';'MLO';'MRO'};
for n = 1:length(changroups)
   subplot(5,2,n)
   plot(erpComm.time,mean(erpICorr.avg(strncmp(erpComm.label,changroups{n},3),:) ,1))
   hold on
   plot(erpComm.time,mean(erpComm.avg(strncmp(erpComm.label,changroups{n},3),:) ,1))
   plot(erpComm.time,mean(erpDiff.avg(strncmp(erpComm.label,changroups{n},3),:) ,1),'k')
   title(changroups{n}); grid on; %ylim([-1 1]*1.5e-13)
end

cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.interactive = 'no';
cfg.showoutline = 'yes';
figure
ft_multiplotER(cfg, erpDiff)
set(gcf,'color','w','position',[ 10   10  1055   849])


% m = permute([erpCCorr.trial;erpICorr.trial;erpComm.trial],[2,1,3]);
% [coeff,score] = pca(reshape(m,[size(m,1),size(m,2)*size(m,3)])');
[coeff,score] = pca([erpCCorr.avg,erpICorr.avg,erpComm.avg]');


figure; set(gcf,'color','w','position',[  683     4   611   876])
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.parameter = 'avg';
cfg.interpolatenan = 'no';
%     cfg.zlim =clim;
cfg.comment    = 'no';
for cc= 1:5
    pcacomp = erpCCorr;
    pcacomp.avg = repmat(coeff(:,cc),[1,length(pcacomp.time)]);

    subplot(5,2,(cc-1)*2 + 1)
    ft_topoplotER(cfg, pcacomp)
    title(sprintf('PC %d',cc))
    
    subplot(5,2,(cc-1)*2 + 2)
    hold on
    if strcmp(stimname,'resp' )
        plot(erpCCorr.time, mean(erpICorr.avg .* coeff(:,cc),1),'linewidth',2)       
        plot(erpComm.time, mean(erpComm.avg .* coeff(:,cc),1),'linewidth',2)
    else
        plot(erpComm.time, mean(erpCCorr.avg .* coeff(:,cc),1),'linewidth',2)
        plot(erpCCorr.time, mean(erpICorr.avg .* coeff(:,cc),1),'linewidth',2)
    end
%     plot(erpComm.time, mean(erpCCorr.avg .* coeff(:,cc),1))
    plot(erpComm.time, mean(erpDiff.avg .* coeff(:,cc),1),'k')
    grid on
end

% if strcmp(stimname,'resp' )
%     saveas(gcf,sprintf('%s/Incong_resp_sens_pca.jpg',resultsfolder))
% else
%     saveas(gcf,sprintf('%s/Corr_flan_sens_pca.jpg',resultsfolder))
% end


%% MNE

gridl =mniLeadfields_singleshell(filenames{1},processingfolder,20484,mri); % calculate leadfields on MNI grid

load([processingfolder,'/headmodel_singleshell.mat']);

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-inf 0]; %it will calculate the covariance matrix
                                   % on the timepoints that are
                                   % before the zero-time point in the trials                             
tlckFC = ft_timelockanalysis(cfg, dataIcomm);
tlckFIC = ft_timelockanalysis(cfg, dataIcorr);

cfg               = [];
cfg.method        = 'mne';
cfg.grid          = gridl;
cfg.headmodel     = vol;
cfg.mne.prewhiten = 'yes';
cfg.mne.lambda    = 3;
cfg.mne.scalesourcecov = 'yes';
sourceFC          = ft_sourceanalysis(cfg,tlckFC);
sourceFIC         = ft_sourceanalysis(cfg, tlckFIC);


figure; set(gcf,'position', [ 46         767        1795         437]) 
clf; t = [0.037, 0.135, 0.15 ];
for jj = 1:3
[~,iit] = min(abs(sourceFC.time - t(jj)));
m=sourceFC.avg.pow(:, iit); % plotting the result at the 450th time-point that is
                         % 500 ms after the zero time-point
subplot(1,3,jj)
ft_plot_mesh(sourceFC, 'vertexcolor', m); caxis([0 1]*4e-22)
view([150 20]); h = light; set(h, 'position', [0 1 0.2]); lighting gouraud; material dull
end

% 
% 
cfg = [];
cfg.projectmom = 'yes';
sdFC  = ft_sourcedescriptives(cfg,sourceFC);
sdFIC = ft_sourcedescriptives(cfg, sourceFIC);

sdDIFF         = sdFC;
sdDIFF.avg.pow = sdFC.avg.pow - sdFIC.avg.pow;

for n = 1:size(sdDIFF.avg.pow,1)
   sdDIFF.avg.pow(n,:) = smooth(sdDIFF.avg.pow(n,:),10) ;
end

cfg = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg,sdDIFF);

%% Co-register MRI

gridres = 5; % 5mm grid
mri_name = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-01/anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_file =  [datapath(1:end-4),'anat/',fids_name];
mri = fids2ctf(mri_name,fids_file,0);

gridl =mniLeadfields_multiSpheres(filenames{1},processingfolder,gridres,mri); % calculate leadfields on MNI grid

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
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
mu  =4;
Cr = C + mu*noiseC; % old normalization
%     Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value


L = gridl.leadfield(gridl.inside);


W = cell(size(L));
for ii = 1:length(L)
 
    lf = L{ii}; % Unit 1Am
    
    % G O'Neill method, equivalent to fieldtrip
    [v,d] = svd(lf'/Cr*lf);
    d = diag(d);
    jj = 2;
    
    lfo = lf*v(:,jj); % Lead field with selected orientation
    
    % no depth correction as we later divide by noise
    w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
    
%     w = Cr\lfo / (lfo'/Cr*lfo) ; % weights
    W{ii} = w;
        
    if mod(ii,300) == 0
        clc
        fprintf('SAM running %.1f\n',...
            ii/length(L)*100)
    end
    
end
%% Beamfomer for PCA


for cc = 1:5
    P = cell(size(L));
    for ii = 1:length(L)
        
        w = W{ii};
        P{ii} = (w'*coeff(:,cc));
  
    end
    
    P  = cell2mat(P)';
    
    crang = [];
    T = zeros(gridl.dim);
    T(gridl.inside) = abs(P);
    sourceant =[];
    sourceant.pow = T;
    sourceant.dim = gridl.dim;
    sourceant.inside = gridl.inside;
    sourceant.pos = gridl.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    % sourceout_Int.coordsys = 'mni';
    
    % crang = [thresh max(sourceant.pow)];
    cfg = [];
    cfg.method        = 'ortho'; %'ortho'
    cfg.location   = 'max';
    
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
%     cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
    
    ft_sourceplot(cfg, sourceout_Int);
    set(gcf,'position',[ 324   270   659   551],'color','w','name',['PC ',num2str(cc)])
    
    
    if strcmp(stimname,'resp' )
        saveas(gcf,sprintf('%s/Incong_%s_beamformed_pca%d.jpg',resultsfolder,stimname,cc))
        
    else
        saveas(gcf,sprintf('%s/Corr_%s_beamformed_pca%d.jpg',resultsfolder,stimname,cc))
    end
    
end



%% Beamfomer


tt= 0.118; 
[~,iit] = min(abs(erpComm.time -tt));
erpIcorrav =  erpICorr.avg(:,iit );
erpcommav=  erpComm.avg(:,iit);
erpCcorrav =  erpCCorr.avg(:,iit);

% erpCorr = ft_timelockanalysis(cfg, datacorr);
% tt= 0.3; 
% [~,iit] = min(abs(erpComm.time -tt));
% erpIcorrav =  erpICorr.avg(:,iit );
% erpcommav=  erpComm.avg(:,iit);
% erpCcorrav =  erpCorr.avg(:,iit);


% timew = [0.25 1];
% % dataA = define_trials([respRIcomm;respLIcomm],data,timeall,timew,1);
% % dataC = define_trials([respRIcomm;respLIcomm],data,timeall,timew-1,1);
% dataA = define_trials([respRIcorr;respLIcorr],data,timeall,timew,1);
% dataC = define_trials([respRIcorr;respLIcorr],data,timeall,timew-1,1);
% 
% Ca = cov(cell2mat(dataA.trial)');
% Cc = cov(cell2mat(dataC.trial)');
% Cc = cov(cell2mat(dataC.trial(randperm(length(dataC.trial),length(dataA.trial))))');




% dataIcorr = define_trials([respRIcorr;respLIcorr],data,timeall,[0.5 1],1);
% dataIcomm = define_trials([respRIcomm;respLIcomm],data,timeall,[0.5 1],1);
% Ca = cov(cell2mat(dataIcomm.trial)');
% % Cc = cov(cell2mat(dataC.trial)');
% Cc = cov(cell2mat(dataIcorr.trial(1:length(dataIcomm.trial)))');

% timew = [0.5 1];
% dataA = define_trials([respRIcorr;respLIcorr],data,timeall,timew,1);
% 
% % dataA = define_trials([respRIcomm;respLIcomm],data,timeall,timew,1);
% dataC = define_trials([respRIcorr;respLIcorr],data,timeall,timew-0.5,1);
% Ca = cov(cell2mat(dataA.trial)');
% % Cc = cov(cell2mat(dataC.trial)');
% Cc = cov(cell2mat(dataC.trial(1:length(dataA.trial)))');

% timew = [0.5 0.2];
% dataA = define_trials([flanRIcomm;flanLIcomm],data,timeall,timew,1);
% dataC = define_trials([flanRIcorr;flanLIcorr],data,timeall,timew,1);
% Ca = cov(cell2mat(dataA.trial)');
% % Cc = cov(cell2mat(dataC.trial)');
% Cc = cov(cell2mat(dataC.trial(1:length(dataA.trial)))');

% timew = [-0.05 0.1];
% dataA = define_trials([buttonresp.right;buttonresp.left],data,timeall,timew,1);
% dataC = define_trials([buttonresp.right;buttonresp.left],data,timeall,timew - 0.2,1);
% Ca = cov(cell2mat(dataA.trial)');
% Cc = cov(cell2mat(dataC.trial)');

% timew = [ 0.12 0.15];
% erpIcorrav =  mean(erpICorr.avg(:,erpComm.time>timew(1) & erpComm.time<timew(2)),2);
% erpcommav=  mean(erpComm.avg(:,erpComm.time>timew(1) & erpComm.time<timew(2)),2);
% erpCcorrav =  mean(erpCCorr.avg(:,erpComm.time>timew(1) & erpComm.time<timew(2)),2);


P = cell(size(L));
Pc = cell(size(L));


for ii = 1:length(L)
    
    w = W{ii};
    
    if strcmp(stimname, 'flan')
        P{ii} = (w'*erpIcorrav); % needs depth corrected weights
        Pc{ii} = (w'*erpCcorrav);
    else
        Pc{ii} = (w'*erpIcorrav); % needs depth corrected weights
        P{ii} = (w'*erpcommav);
    end
    
%         P{ii} = ((w'*Ca*w) - (w'*Cc*w) ) / (2*w'*noiseC*w);
%         P{ii} = (w'*Ca*w)  / (w'*noiseC*w);
    %     Pc{ii} = (w'*Cc*w)  / (w'*noiseC*w);
    
    
end

P  = cell2mat(P)';
Pc  = cell2mat(Pc)';

%%
mriplot = 'brain';
plotopt = 'ortho';  % 'slice'; %'ortho'


crang = [];

T = zeros(gridl.dim);
T(gridl.inside) = abs(P) - abs(Pc) ;
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

% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = plotopt; %'ortho'
cfg.location   = 'max';

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

% set(gcf,'name',sprintf('Corr ERF peak at %dms ',round(tt*1e3)))
% saveas(gcf,sprintf('%s/Comm_ERF_%dms_%s.jpg',resultsfolder,round(tt*1e3),stimname))

% title(sprintf('Incong Comm %d,%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm_%d-%dms_resp.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000))
% 
% saveas(gcf,sprintf('%s/Incong_Corr-base_%d-%dms_13-30Hz_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000,stimname))



%% For ERPs
cc = [];
ccd  =[-3 3];
crang = [cc]*1e-14;

T = zeros(gridl.dim);
T(gridl.inside) = abs(P);
sourceant =[];
sourceant.pow = T;
sourceant.dim = gridl.dim;
sourceant.inside = gridl.inside;
sourceant.pos = gridl.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
% sourceout_Int.coordsys = 'mni';

% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = 'slice'; %'ortho'
cfg.location   = 'max';

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);
set(gcf,'position',[704   117   894   675],'color','w')

title(sprintf('Incong Comm %d-%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm_%d-%dms_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000,stimname))


T = zeros(gridl.dim);
T(gridl.inside) = abs(Pc);%-abs(Pc);
sourceant =[];
sourceant.pow = T;
sourceant.dim = gridl.dim;
sourceant.inside = gridl.inside;
sourceant.pos = gridl.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
% sourceout_Int.coordsys = 'mni';

% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = 'slice'; %'ortho'
cfg.location   = 'max';

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);
set(gcf,'position',[704   117   894   675],'color','w')
% title(['Incong Comm - Corr 200ms after response'])
title(sprintf('Incong Corr %d-%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absCorr_%d-%dms_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000, stimname))


crang = [cc]*1e-14;
T = zeros(gridl.dim);
T(gridl.inside) = abs(P-Pc);%-abs(Pc);
sourceant =[];
sourceant.pow = T;
sourceant.dim = gridl.dim;
sourceant.inside = gridl.inside;
sourceant.pos = gridl.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
% sourceout_Int.coordsys = 'mni';

% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = 'slice'; %'ortho'
cfg.location   = 'max';

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);
set(gcf,'position',[704   117   894   675],'color','w')
title(sprintf('abs(Incong Comm - Corr) %d-%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm-Corr_%d-%dms_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000, stimname))




crang = [ccd]*1e-14;
T = zeros(gridl.dim);
T(gridl.inside) = abs(P)-abs(Pc);
sourceant =[];
sourceant.pow = T;
sourceant.dim = gridl.dim;
sourceant.inside = gridl.inside;
sourceant.pos = gridl.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
% sourceout_Int.coordsys = 'mni';

% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = 'slice'; %'ortho'
cfg.location   = 'max';

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);
set(gcf,'position',[704   117   894   675],'color','w')
title(sprintf('abs(Incong Comm) - abs(Incong Corr) %d-%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm-absCorr_%d-%dms_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000, stimname))



end

