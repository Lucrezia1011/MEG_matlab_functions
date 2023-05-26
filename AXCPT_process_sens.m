%function AXCPT_process_sens(filenames,datapath,processingfolder,bvfiles,fids_name)
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
%% Standard pre-processing

cd(datapath)

highpass = 0.5;
lowpass = 120;
icaopt = 1;
plotopt = 0;
realignopt = 1;

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
        
        if realignopt == 1
            if ~exist(mri_name,'file')
                mri_name = [mri_name,'.gz'];
            end
            fids_file =  [datapath(1:end-4),'anat/',fids_name];
            mri = fids2ctf(mri_name,fids_file,0);

            mneres =5124; % 5124, 8196, 20484
            gridl =mniLeadfields_singleshell(filename,processingfolder,mneres,mri); % calculate leadfields on MNI grid
            load([processingfolder,'/headmodel_singleshell.mat']);
        end
    else
        
        if any(strcmp(data.hdr.label,'UADC009'))
            eyeall = cat(2,eyeall,eyelink);
        end
        if realignopt == 1
         % realign sensors to first recording
            cfg = [];
            cfg.template          = dataall.grad;
            cfg.headmodel  =  vol;
            cfg.inwardshift = 1; 
            [data] = ft_megrealign(cfg, data);
         end
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

close all
%% Filters
filt_order = []; % default
data_filt = ft_preproc_lowpassfilter(dataall.trial{1}, dataall.fsample,50,filt_order,'but');
% data_filt = ft_preproc_bandpassfilter(dataall.trial{1}, dataall.fsample,[13 30],filt_order,'but');
data = dataall;
data.trial{1} = data_filt;
clear data_filt


%% Pupil
pupil = eyeall(3,:);

[pks,locs,w,p] = findpeaks(-pupil,'MinPeakDistance',6,'MinPeakProminence',0.6);
% figure
% plot(1:length(pupil),pupil,(locs),-pks,'or')
% xlabel('time')
% ylabel('pupil')
% title('blinks')
% axis tight
% legend('Data','peaks','Location','NorthWest')

pupiln = pupil;
buff = 30;
for n = 1:length(locs)
    stepf = 1;
    df = pupil(locs(n)+stepf) - pupil(locs(n));
    
    while df < p(n)*0.8
        stepf = stepf + 1;
        df = pupil(locs(n)+stepf) - pupil(locs(n));
    end
    
    stepb = 1;
    db = pupil(locs(n)-stepb) - pupil(locs(n));
    
    while db < p(n)*0.8
        stepb = stepb + 1;
        db = pupil(locs(n)-stepb) - pupil(locs(n));
    end
    
    pupiln( (locs(n)-stepb-buff): (locs(n)+stepf+buff) ) = NaN;
%     pupil( (locs(n)-ceil(buff*p(n))): (locs(n)+ceil(w(n)+buff*p(n))) ) = NaN;
    
end
% plot(pupil)
% hold on
% plot(pupiln)
% 
pupili = fillmissing(pupiln,'pchip');
% plot(pupili)


data.label(end+1:end+3) = {'eyeX';'eyeY';'pupil'};
data.trial{1}(end+1:end+3,:) = [eyeall(1:2,:);pupili];



%% Define Trials

twind = [-1 2];
cue = 'A';
% buttonpress to cue
cuesampA = trigall.sample(strcmp(trigall.type, cue) ) ;

buttonsamp = match_responses(cuesampA, buttonresp, 'left', data.fsample);
[datacueAresp,ttdel] = define_trials(buttonsamp,data,timeall,twind,1);


cue = 'B';
% buttonpress to cue
cuesampB = trigall.sample(strcmp(trigall.type, cue) ) ;


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

probe = 'BY';
%buttonpress to probe
probesampBY = trigall.sample(strcmp(trigall.type, probe)) ;


buttonsampBXcomm = match_responses(probesampBX, buttonresp, 'right', data.fsample);
buttonsampBXcorr = match_responses(probesampBX, buttonresp, 'left', data.fsample);

[datacueA,ttdel] = define_trials(cuesampA,data,timeall,twind,1);
[datacueB,ttdel] = define_trials(cuesampB,data,timeall,twind,1);


[dataprobeAX,ttdel] = define_trials(probesampAX,data,timeall,twind,1);
[dataprobeBX,ttdel] = define_trials(probesampBX,data,timeall,twind,1);
[dataprobeBY,ttdel] = define_trials(probesampBY,data,timeall,twind,1);
[dataprobeAY,ttdel] = define_trials(probesampAY,data,timeall,twind,1);

[dataproberespComm,ttdel] = define_trials([buttonsampAYcomm;buttonsampAXcomm],data,timeall,twind,1);
[dataproberespAYComm,ttdel] = define_trials([buttonsampAYcomm],data,timeall,twind,1);

[dataproberespAXCorr,ttdel] = define_trials(buttonsampAXcorr,data,timeall,twind,1);
[dataproberespAYCorr,ttdel] = define_trials(buttonsampAYcorr,data,timeall,twind,1);
[dataprobeArespCorr,~] =  define_trials([buttonsampAXcorr;buttonsampAYcorr],data,timeall,twind,1);

cfg = [];
cfg.demean        = 'yes';
cfg.baselinewindow = [-inf, 0];

if strcmp(cfg.demean, 'yes')
datacueA = ft_preprocessing(cfg,datacueA);
datacueB = ft_preprocessing(cfg,datacueB);

dataprobeAX = ft_preprocessing(cfg,dataprobeAX);
dataprobeBX = ft_preprocessing(cfg,dataprobeBX);
dataprobeAY = ft_preprocessing(cfg,dataprobeAY);
dataproberespComm = ft_preprocessing(cfg,dataproberespComm);
dataproberespAYComm = ft_preprocessing(cfg,dataproberespAYComm);

dataproberespAXCorr = ft_preprocessing(cfg,dataproberespAXCorr);
dataproberespAYCorr = ft_preprocessing(cfg,dataproberespAYCorr);
end



erpComm = ft_timelockanalysis([], dataproberespComm);
erpCommAY = ft_timelockanalysis([], dataproberespAYComm);

erpCorrAX = ft_timelockanalysis([], dataproberespAXCorr);
erpCorrAY = ft_timelockanalysis([], dataproberespAYCorr);

erpProbeAX = ft_timelockanalysis([], dataprobeAX);
erpProbeAY = ft_timelockanalysis([], dataprobeAY);

erpProbeBX = ft_timelockanalysis([], dataprobeBX);
erpProbeBY = ft_timelockanalysis([], dataprobeBY);

erpCueA = ft_timelockanalysis([], datacueA);
erpCueB = ft_timelockanalysis([], datacueB);

%%
figure(1); set(gcf,'color','w')
clf
changroups = {'MLF';'MRF';'MLC';'MRC';'MLT';'MRT';'MLP';'MRP';'MLO';'MRO';'eye';'pup'};
for n = 1:length(changroups)
   subplot(6,2,n)
   hold on
   plot(erpProbeAX.time,mean(erpProbeAX.avg(strncmp(erpProbeAX.label,changroups{n},3),:) ,1),'linewidth',2)
   plot(erpProbeAY.time,mean(erpProbeAY.avg(strncmp(erpProbeAY.label,changroups{n},3),:) ,1),'linewidth',2)
%    plot(erpComm.time,mean(erpComm.avg(strncmp(erpComm.label,changroups{n},3),:) ,1))
   xlim([-0.5 1])
   title(changroups{n}); grid on; %ylim([-1 1]*1.5e-13)
end
legend('AX','AY')

figure(2);set(gcf,'color','w')
clf
for n = 1:length(changroups)
   subplot(6,2,n)
   hold on
   plot(erpComm.time,mean(erpCorrAX.avg(strncmp(erpCorrAX.label,changroups{n},3),:) ,1),'linewidth',2)
   plot(erpComm.time,mean(erpCorrAY.avg(strncmp(erpCorrAY.label,changroups{n},3),:) ,1),'linewidth',2)
%    plot(erpComm.time,mean(erpComm.avg(strncmp(erpComm.label,changroups{n},3),:) ,1))
   xlim([-0.5 1])
   title(changroups{n}); grid on; %ylim([-1 1]*1.5e-13)
end
legend('AX corr','AY corr')


figure(3);set(gcf,'color','w')
clf
for n = 1:length(changroups)
   subplot(6,2,n)
   hold on
   plot(erpCueA.time,mean(erpCueA.avg(strncmp(erpCorrAX.label,changroups{n},3),:) ,1),'linewidth',2)
   plot(erpCueB.time,mean(erpCueB.avg(strncmp(erpCorrAY.label,changroups{n},3),:) ,1),'linewidth',2)
%    plot(erpComm.time,mean(erpComm.avg(strncmp(erpComm.label,changroups{n},3),:) ,1))
   xlim([-0.5 1])
   title(changroups{n}); grid on; %ylim([-1 1]*1.5e-13)
end
legend('A','B')

%% ERP topoplot over time
figure; set(gcf,'color','w','position',[  32     481   1771    399])
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.parameter = 'avg';
cfg.interpolatenan = 'no';
cfg.colorbar = 'SouthOutside'  ;
cfg.zlim =[-1 1]*5e-14;
cfg.comment    = 'no';
for cc= 1:6
subplot(1,6,cc)
cfg.xlim = 0.025 + [(cc-1),cc]*0.05; % time 
ft_topoplotER(cfg, erpProbeAY)
title(sprintf('%d-%dms',round(cfg.xlim(1)*1e3),round(cfg.xlim(2)*1e3)))
end
%% 
erp1 = erpCueA;
erp2 = erpCueB;

% erp1 = erpProbeAX;
% erp2 = erpProbeAY;

mchans = strncmp(data.label,'M',1);
[coeff,score] = pca([erp1.avg(mchans,:),erp2.avg(mchans,:)]');

figure; set(gcf,'color','w','position',[  683     4   611   876])
clf
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.parameter = 'avg';
cfg.interpolatenan = 'no';
co = get(gca,'colororder');
%     cfg.zlim =clim;
cfg.comment    = 'no';
for cc= 1:5
    pcacomp = erp1;
    pcacomp.var(~mchans,:) = [];
    pcacomp.dof(~mchans,:) = [];
    pcacomp.label(~mchans) = [];
    pcacomp.avg = repmat(coeff(:,cc),[1,length(pcacomp.time)]);
    
    subplot(5,3,(cc-1)*3 + 1)
    ft_topoplotER(cfg, pcacomp)
    title(sprintf('PC %d',cc))
    
    subplot(5,3,(cc-1)*3 + (2:3))
    hold on
    
    plot(erp1.time, mean(erp1.avg(mchans,:) .* coeff(:,cc),1),'linewidth',2)
    plot(erp1.time, mean(erp2.avg(mchans,:) .* coeff(:,cc),1),'linewidth',2)
 
    xlim([-0.5 1])
%     plot(erp1.time, mean(erp1.avg(mchans,:) .* coeff(:,cc),1),'linewidth',2,'color',co(1,:))
%     plot(erp1.time, mean(erp2.avg(mchans,:) .* coeff(:,cc),1),'--','linewidth',2,'color',co(1,:))
%     plot(erp1.time, mean(erp3.avg(mchans,:) .* coeff(:,cc),1),'linewidth',2,'color',co(2,:))
%     plot(erp1.time, mean(erp4.avg(mchans,:) .* coeff(:,cc),1),'--','linewidth',2,'color',co(2,:))
%     y = ylim;
%     plot([0.3 0.3],y,'k--')
    %     plot(erpComm.time, mean(erpDiff.avg .* coeff(:,cc),1),'k')
    grid on
end
% legend('AX','AY','location','best')
% saveas(gcf,sprintf('%s/Probe_ERF_sens_05-50Hz.jpg',resultsfolder))

legend('A','B','location','best')
saveas(gcf,sprintf('%s/Cue_ERF_sens_05-50Hz.jpg',resultsfolder))



%% TFS
changroups = {'MRT','MLT','MRF','MLF','MRO','MLO','MRC','MLC','MRP','MLP'};

cfg              = [];
cfg.output       = 'pow';
cfg.channel =   'MEG';%datacue.label(strncmp(datacue.label,changroups{c},3));
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:2:50;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.toi          = -1:0.05:2;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
TFRhannAX = ft_freqanalysis(cfg, dataprobeAX );
TFRhannAY = ft_freqanalysis(cfg, dataprobeAY);
TFRhannA = ft_freqanalysis(cfg, datacueA );
TFRhannB = ft_freqanalysis(cfg, datacueB);

%%
figure; set(gcf,'color','w','position', [1008     5     615     870])
cfg                 = [];
cfg.baseline        = [-1 -0.3];
cfg.baselinetype    = 'relchange';
cfg.zlim            = [-0.8 0.8]; %'maxabs';
% cfg.ylim            = [0 50];
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,TFRhannAX)
    title(changroups{n});
end
saveas(gcf,sprintf('%s/ProbeAX_TFS_sens.jpg',resultsfolder))

figure; set(gcf,'color','w','position', [1008     5     615     870])
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,TFRhannAY)
    title(changroups{n});
end
saveas(gcf,sprintf('%s/ProbeAY_TFS_sens.jpg',resultsfolder))


figure; set(gcf,'color','w','position', [1008     5     615     870])
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,TFRhannA)
    title(changroups{n});
end
saveas(gcf,sprintf('%s/CueA_TFS_sens.jpg',resultsfolder))


figure; set(gcf,'color','w','position', [1008     5     615     870])
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,TFRhannB)
    title(changroups{n});
end
saveas(gcf,sprintf('%s/CueB_TFS_sens.jpg',resultsfolder))

TFRhann = TFRhannA;
TFRhann.powspctrm = TFRhannAX.powspctrm - TFRhannAY.powspctrm;

figure; set(gcf,'color','w','position', [1008     5     615     870])
cfg                 = [];
% cfg.baseline        = [-1 -0.3];
% cfg.baselinetype    = 'absolute';%'relchange';
cfg.zlim            = [-1 1]*1e-27; %'maxabs';
% cfg.ylim            = [0 50];
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,TFRhann)
    title(changroups{n});
end

saveas(gcf,sprintf('%s/ProbeAX-AY_TFS_sens.jpg',resultsfolder))


TFRhann.powspctrm = TFRhannA.powspctrm - TFRhannB.powspctrm;
figure; set(gcf,'color','w','position', [1008     5     615     870])
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,TFRhann)
    title(changroups{n});
end

saveas(gcf,sprintf('%s/CueA-B_TFS_sens.jpg',resultsfolder))


