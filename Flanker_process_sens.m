% function Flanker_process_sens(filenames,datapath,processingfolder,bvfiles,fids_name)

addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

sub = '24580';

datapath = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/meg/'];

processingfolder = ['/data/liuzzil2/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
if ~exist(processingfolder,'dir')
    mkdir(processingfolder)
end

 

filenames = cell(1,2);
for ii = 1:2
    filenames{ii} = ['sub-',sub,'_task-flanker_run-',num2str(ii),'_meg.ds'];
end
mri_name = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/anat/sub-',sub,'_acq-mprage_T1w.nii'];
fids_name = ['sub-',sub,'_fiducials.tag'];

%% Standard pre-processing

cd(datapath)
 
highpass = 0.5;
lowpass = 120;
icaopt = 1;
plotopt = 0;

realignopt = 1;

for ii = 1:length(filenames)
    filename = filenames{ii};
%     sub = filename(5:9);
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
    trig_sample.sample(trig_sample.sample == 0) = [];
    if ~isempty(buttonpress)
        buttonpress.left = samples_T(buttonpress.UADC006);
        buttonpress.left(buttonpress.left == 0) = [];
        buttonpress.right = samples_T(buttonpress.UADC007);
        buttonpress.right(buttonpress.right == 0) = [];
    else
        % Read behavioral file
        M = readtable(['../behavioral/',bvfiles{ii}]) ;
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
        if realignopt == 1
            if ~exist(mri_name,'file')
                mri_name = [mri_name,'.gz'];
            end
            fids_file =  [datapath(1:end-4),'anat/',fids_name];
            mri = fids2ctf(mri_name,fids_file,0);

            mneres =20484; % 5124, 8196, 20484
            gridl =mniLeadfields_singleshell(filename,processingfolder,mneres,mri); % calculate leadfields on MNI grid
            load([processingfolder,'/headmodel_singleshell.mat']);
        end
        
    else
        
        if any(strcmp(data.hdr.label,'UADC009'))
            eyeall = cat(2,eyeall,eyelink);
        end
        % realign sensors to first recording
        if realignopt == 1
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



%% ICA of concatenated data
resultsfolder = [processingfolder,'flanker/'];

if icaopt == 0

    rmpath(('~/fieldtrip-20190812/fieldtrip_private'))
    if exist([resultsfolder,'/ICA_artifacts.mat'],'file')
        load([resultsfolder,'/ICA_artifacts.mat']);
    else
        if ~exist(resultsfolder,'dir')
            mkdir(resultsfolder)
        end
        dataica = dataall;
        if highpass < 1
            dataica.trial{1} = ft_preproc_highpassfilter(dataall.trial{1}, dataall.fsample, 1, 4, [], [], []);     
        end

        cfg =[];
        cfg.method = 'pca';
        comp_pca = ft_componentanalysis(cfg, dataica);
        score = comp_pca.trial{1}';
        compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
        icomp = nnz(compvar95) ;
        clc
        fprintf('%d components for 95perc. of data variance\n',icomp)

        if icomp>30
            disp('Reducing ICA components to 30')
            icomp = 30;
        end
        cfg =[];
        cfg.method = 'fastica';
        cfg.fastica.numOfIC = icomp;
        cfg.fastica.maxNumIterations = 100;
        %         cfg.icasso.mode
        comp = ft_componentanalysis(cfg, dataica);
        icomp = length(comp.label);

        %         figure
        %         cfg           = [];
        %         cfg.component = [1:icomp];       % specify the component(s) that should be plotted
        %         cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
        %         cfg.comment   = 'no';
        %         ft_topoplotIC(cfg, comp)


        cfg          = [];
        cfg.channel  = [1:5]; % components to be plotted
        cfg.viewmode = 'component';
        cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
        cfg.continuous   = 'yes' ;
        cfg.blocksize    = 10;
        ft_databrowser(cfg, comp)

        if exist('eyeall','var')
            figure; % delete bad data segments
            plot(abs(corr(eyeall',comp.trial{1}'))','*')
            grid on
            xlabel('ICA component')
            ylabel('Correlation with Eye movements')
        end
        icadel = input('ICA component to eliminate (input as [''01'';''02'']): ');

        cfg = [];
        cfg.channel = cell(size(icadel,1),1);
        for ii = 1:size(icadel,1)
            cfg.channel{ii}  = ['fastica0',icadel(ii,:)];
        end

        [comps] = ft_selectdata(cfg, comp);
        save([resultsfolder,'/ICA_artifacts'],'comps')
        close all
    end

    cfg           = [];
    cfg.component = 1:length(comps.label);
    dataall          = ft_rejectcomponent(cfg, comps,dataall);
    
end

%% 
filt_order = []; % default
data_filt = ft_preproc_lowpassfilter(dataall.trial{1}, dataall.fsample,55,filt_order,'but');
% data_filt = ft_preproc_bandpassfilter(dataall.trial{1}, dataall.fsample,[13 30],filt_order,'but');
data = dataall;
data.trial{1} = data_filt;
clear data_filt


%%

  
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

flan = 'LC';
resp = 'right';
% buttonpress to cue
cuesamp = trigall.sample(strcmp(trigall.type, flan) ) ;
[respLCcomm,flanLCcomm] = resp2trial(resp,cuesamp,buttonresp,data.fsample);

flan = 'RC';
resp = 'left';
% buttonpress to cue
cuesamp = trigall.sample(strcmp(trigall.type, flan) ) ;
[respRCcomm,flanRCcomm] = resp2trial(resp,cuesamp,buttonresp,data.fsample);

%% Plot sens
stimname = 'flan';
if strcmp(stimname,'resp' )
    timewresp = [-1,1];
    bsw = [-1, -0.5];
%     dataIcorr = define_trials([respRIcorr;respLIcorr],data,timeall,timewresp,1);
%     dataIcomm = define_trials([respRIcomm;respLIcomm],data,timeall,timewresp,1);
%     dataCcorr = define_trials([respRCcorr;respLCcorr],data,timeall,timewresp,1);
    dataC = define_trials([respRCcorr;respLCcorr;respRCcomm;respLCcomm],data,timeall,timewresp,1);
    dataI = define_trials([respRIcorr;respLIcorr;respRIcomm;respLIcomm],data,timeall,timewresp,1);
else
    timewresp = [-0.5,1.5];
    bsw = [-0.5, -0];
%     dataIcorr = define_trials([flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
%     dataIcomm = define_trials([flanRIcomm;flanLIcomm],data,timeall,timewresp,1);
%     dataCcorr = define_trials([flanRCcorr;flanLCcorr],data,timeall,timewresp,1);
    dataC = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],data,timeall,timewresp,1);
    dataI = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],data,timeall,timewresp,1);
%     datacorr = define_trials([flanRCcorr;flanLCcorr;flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
end


cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = bsw;
% dataIcomm = ft_preprocessing(cfg,dataIcomm);
% dataIcorr = ft_preprocessing(cfg,dataIcorr);
% dataCcorr = ft_preprocessing(cfg,dataCcorr);
dataI = ft_preprocessing(cfg,dataI);
dataC = ft_preprocessing(cfg,dataC);

cfg= [];
cfg.keeptrials         = 'no';
% erpComm = ft_timelockanalysis(cfg, dataIcomm);
% erpICorr = ft_timelockanalysis(cfg, dataIcorr);
% erpCCorr = ft_timelockanalysis(cfg, dataCcorr);

erpI = ft_timelockanalysis(cfg, dataI);
erpC = ft_timelockanalysis(cfg, dataC);

% if strcmp(stimname,'resp' )
%     erpDiff = erpComm;
%     erpDiff.avg = erpComm.avg - erpICorr.avg;
% else
%     erpDiff = erpICorr;
%     erpDiff.avg = erpICorr.avg - erpCCorr.avg;
% end
erpDiff = erpI;
erpDiff.avg = erpI.avg - erpC.avg;

figure; set(gcf,'color','w','position',[647    67   667   767])
clf
changroups = {'MLF';'MRF';'MLC';'MRC';'MLT';'MRT';'MLP';'MRP';'MLO';'MRO'};
for n = 1:length(changroups)
   subplot(5,2,n)
   hold all
%    if strcmp(stimname,'resp' )
%         plot(erpCCorr.time, mean(erpICorr.avg(strncmp(erpComm.label,changroups{n},3),:) ,1),'linewidth',2)       
%         plot(erpComm.time, mean(erpComm.avg(strncmp(erpComm.label,changroups{n},3),:) ,1),'linewidth',2)
%     else
%         plot(erpComm.time, mean(erpCCorr.avg(strncmp(erpComm.label,changroups{n},3),:) ,1),'linewidth',2)
%         plot(erpCCorr.time, mean(erpICorr.avg(strncmp(erpComm.label,changroups{n},3),:),1),'linewidth',2)
%    end
   plot(erpC.time, mean(erpC.avg(strncmp(erpComm.label,changroups{n},3),:) ,1),'linewidth',2)
   plot(erpC.time, mean(erpI.avg(strncmp(erpComm.label,changroups{n},3),:),1),'linewidth',2)
    
   plot(erpC.time,mean(erpDiff.avg(strncmp(erpComm.label,changroups{n},3),:) ,1),'k')
   title(changroups{n}); grid on; %ylim([-1 1]*1.5e-13)
%    xlim([-0.5 1])
end

legend('cong','incg','diff','location','best')
%%
% if strcmp(stimname,'resp' )
% %     legend('correct','commission','diff','location','best')
%     saveas(gcf,sprintf('%s/Incong_resp_sens_05hzhp.jpg',resultsfolder))
% else
% %     legend('congruent','incongruent','diff','location','best')
%     saveas(gcf,sprintf('%s/Corr_flan_sens_05hzhp.jpg',resultsfolder))
% end

% cfg = [];
% cfg.layout = 'CTF275_helmet.mat';
% cfg.interactive = 'no';
% cfg.showoutline = 'yes';
% figure
% ft_multiplotER(cfg, erpDiff)
% set(gcf,'color','w','position',[ 10   10  1055   849])


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

if strcmp(stimname,'resp' )
    saveas(gcf,sprintf('%s/Incong_resp_sens_pca_01hzhp.jpg',resultsfolder))
else
    saveas(gcf,sprintf('%s/Corr_flan_sens_pca_01hzhp.jpg',resultsfolder))
end

%% TFS
close all

timewresp = [-1 2];   

stimname = 'resp';
if strcmp(stimname,'resp' )
%     dataIcorr = define_trials([respRIcorr;respLIcorr],data,timeall,timewresp,1);
%     dataIcomm = define_trials([respRIcomm;respLIcomm],data,timeall,timewresp,1);
%     dataCcorr = define_trials([respRCcorr;respLCcorr],data,timeall,timewresp,1);
    dataC = define_trials([respRCcorr;respLCcorr;respRCcomm;respLCcomm],data,timeall,timewresp,1);
    dataI = define_trials([respRIcorr;respLIcorr;respRIcomm;respLIcomm],data,timeall,timewresp,1);

    
else
%     dataIcorr = define_trials([flanRIcorr;flanLIcorr],data,timeall,timewresp,1);
%     dataIcomm = define_trials([flanRIcomm;flanLIcomm],data,timeall,timewresp,1);
%     dataCcorr = define_trials([flanRCcorr;flanLCcorr],data,timeall,timewresp,1);
    dataC = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],data,timeall,timewresp,1);
    dataI = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],data,timeall,timewresp,1);

end



cfg                 = [];
cfg.method          = 'mtmconvol';
cfg.taper           = 'hanning';
cfg.channel         = 'MEG';

% set the frequencies of interest
cfg.foi             = 1:2:55;
% set the timepoints of interest: from -0.8 to 1.1 in steps of 100ms
cfg.toi             = -1:0.1:2;

% set the time window for TFR analysis: constant length of 200ms
cfg.t_ftimwin       = 0.2 * ones(length(cfg.foi), 1);

% average over trials
cfg.keeptrials      = 'no';

% pad trials to integer number of seconds, this speeds up the analysis
% and results in a neatly spaced frequency axis
cfg.pad             = 4;
% freqIcorr           = ft_freqanalysis(cfg, dataIcorr);
% freqIcomm           = ft_freqanalysis(cfg, dataIcomm);
% freqCcorr           = ft_freqanalysis(cfg, dataCcorr);
freqI           = ft_freqanalysis(cfg, dataI);
freqC           = ft_freqanalysis(cfg, dataC);

% if strcmp(stimname,'resp' )
%     freqdiff = freqIcorr;
%     freqdiff.powspctrm = freqIcomm.powspctrm - freqIcorr.powspctrm;
% else
%     freqdiff = freqIcorr;
%     freqdiff.powspctrm = freqIcorr.powspctrm - freqCcorr.powspctrm;
% 
% end
freqdiff = freqI;
freqdiff.powspctrm = freqI.powspctrm - freqC.powspctrm;

changroups = {'MLF';'MRF';'MLC';'MRC';'MLT';'MRT';'MLP';'MRP';'MLO';'MRO'};

cfg                 = [];
cfg.baseline        = timewresp;%[-1 -0.3];
cfg.baselinetype    = 'relchange';
cfg.zlim            = [-1 1];
cfg.ylim            = [0 50];

figure; set(gcf,'color','w','position',[135    53   590   866],'name','Congruent')

for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,freqC)
    title(changroups{n});
end
saveas(gcf,sprintf('%s/TFS_Cong_%s_sens.jpg',resultsfolder,stimname))


figure; set(gcf,'color','w','position',[135    53   590   866],'name','Incongruent')
for n = 1:length(changroups)
    subplot(5,2,n)
    cfg.channel = data.label(strncmp(data.label,changroups{n},3));
    ft_singleplotTFR(cfg,freqI)
    title(changroups{n});
end
saveas(gcf,sprintf('%s/TFS_Incong_%s_sens.jpg',resultsfolder,stimname))
% 
% cfg                 = [];
% cfg.zlim            ='maxabs';
% cfg.ylim            = [0 50];
% for n = 1:length(changroups)
%     subplot(5,2,n)
%     cfg.channel = data.label(strncmp(data.label,changroups{n},3));
%     ft_singleplotTFR(cfg,freqdiff)
%     title(changroups{n});
% end
% saveas(gcf,sprintf('%s/TFS_diff_%s_sens.jpg',resultsfolder,stimname))
% 
% 


