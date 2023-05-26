% function Flanker_process(filenames,datapath,processingfolder,bvfiles,fids_name)

addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

sublist = {'24531';'24563';'24590';'24580';'24588';'24626';'24581';'24482';...
    '24640';'24592';'24667';'24678';'24663' };
% 24592 could not extract blink ICA
for ss = length(sublist)
    sub = sublist{ss};
    
    if strcmp(sub,'24531')
        datapath = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-02/meg/'];
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/ses-02/'];
        mri_name = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-01/anat/sub-',sub,'_acq-mprage_T1w.nii'];
    else
        datapath = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sub,'/meg/'];
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
        mri_name = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sub,'/anat/sub-',sub,'_acq-mprage_T1w.nii'];
    end
    
    if ~exist(processingfolder,'dir')
        mkdir(processingfolder)
    end
    
    d = dir(datapath);
    filenames = cell(1,2);
    n = 0;
    for ii = 3:length(d)
        if  contains(d(ii).name,'meg.ds') && contains(d(ii).name,'task-flanker')
            n = n+1;
            filenames{n} = d(ii).name;
        end
        
    end
    
    fids_name = ['sub-',sub,'_fiducials.tag'];
    
    close all
    
    
    cd(datapath)
    
    if ~exist(mri_name,'file')
        mri_name = [mri_name,'.gz'];
    end
    fids_file =  [datapath(1:end-4),'anat/',fids_name];
    if ~exist(fids_file,'file')
        fids_file =  [datapath(1:end-4),'anat/sub-',sub,'_fiducials_flanker.tag'];
    end
    mri = fids2ctf(mri_name,fids_file,0);
    
    
    
    
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
    
    d = dir('/data/EDB/MEG_AXCPT_Flanker/data/emptyroom/');
    emptyroom = []; jj = 0;
    for ii = 3:length(d)
        if contains(d(ii).name, TaskDate(1:8))
            jj = jj + 1;
            emptyroom{jj} = ['/data/EDB/MEG_AXCPT_Flanker/data/emptyroom/',d(ii).name];
        end
    end
    %% Standard pre-processing
    
    highpass = 0.5;
    lowpass = 120;
    icaopt = 1;
    plotopt = 0;
    
    
    % Data header info
    hdr = ft_read_header(emptyroom{1});
    % Get Bad channel names
    fid = fopen([emptyroom{1},'/BadChannels']);
    BadChannels = textscan(fid,'%s');
    fclose(fid);
    
    % get MEG channel names
    channels = hdr.label(strcmp(hdr.chantype,'meggrad'));
    % Delete Bad channels
    chanInd = zeros(size(channels));
    for iiC = 1:length(BadChannels{1})
        chanInd = chanInd | strcmp(channels,BadChannels{1}{iiC});
    end
    channels(find(chanInd)) = [];
    
    noiseC = zeros(length(channels),length(channels),length(emptyroom));
    noiseCtheta = zeros(length(channels),length(channels),length(emptyroom));
    noiseCalpha = zeros(length(channels),length(channels),length(emptyroom));
    
    for ii = 1:length(emptyroom)
        cfg = [];
        cfg.dataset = emptyroom{ii};
        cfg.continuous = 'yes';
        cfg.channel = channels;
        cfg.demean = 'yes';
        cfg.detrend = 'no';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [lowpass, highpass];
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [58 62; 118 122; 178 182]; % With notch filter 60Hz
        
        data_empty = ft_preprocessing(cfg);
        emptyC = zeros(length(channels),length(channels),length(data_empty.trial));
        for t = 1:length(data_empty.trial)
            emptyC = cov(data_empty.trial{t}');
        end
        noiseC(:,:,ii) = mean(emptyC,3);
        
        
        
        cfg.bpfreq = [4, 8];
        data_empty = ft_preprocessing(cfg);
        emptyC = zeros(length(channels),length(channels),length(data_empty.trial));
        for t = 1:length(data_empty.trial)
            emptyC = cov(data_empty.trial{t}');
        end
        noiseCtheta(:,:,ii) = mean(emptyC,3);
        
        
        cfg.bpfreq = [8 13];
        data_empty = ft_preprocessing(cfg);
        emptyC = zeros(length(channels),length(channels),length(data_empty.trial));
        for t = 1:length(data_empty.trial)
            emptyC = cov(data_empty.trial{t}');
        end
        noiseCalpha(:,:,ii) = mean(emptyC,3);
    end
    noiseC = mean(noiseC,3);
    noiseCtheta = mean(noiseCtheta,3);
    noiseCalpha = mean(noiseCalpha,3);
    
    %%
    cd(datapath)
    resultsfolder = [processingfolder,'flanker/'];
    
    
    for ii = 1:length(filenames)
        filename = filenames{ii};
        gridres = 5; % 5mm grid
        mu  =4; % regularization parameter
        if ~exist(sprintf('%s%s/ERN_multiSpheres_%dmm_regmu%d.mat',...
                processingfolder,filename(1:end-3),gridres,mu),'file')
            
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
            trig_sample.type(trig_sample.sample == 0) = [];
            trig_sample.value(trig_sample.sample == 0) = [];
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
            
            
            %     trig_sample = trig_sample;
            %     buttonresp = buttonpress;
            %     dataall = data;
            %     timeall = time;
            %     if any(strcmp(data.hdr.label,'UADC009'))
            %         eyeall = eyelink;
            %     end
            %
            %         mri_name = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-01/anat/sub-',sub,'_acq-mprage_T1w.nii'];
            
            
            %
            %     mneres =20484; % 5124, 8196, 20484
            %     gridl =mniLeadfields_singleshell(filename,processingfolder,mneres,mri); % calculate leadfields on MNI grid
            %     load([processingfolder,'/headmodel_singleshell.mat']);
            %
            
            %%
            
            flan = 'RI';
            resp = 'left';
            % buttonpress to cue
            cuesampr = trig_sample.sample(strcmp(trig_sample.type, flan) ) ;
            [respRIcorr,flanRIcorr] = resp2trial(resp,cuesampr,buttonpress,data.fsample);
            
            resp = 'right';
            [respRIcomm,flanRIcomm] = resp2trial(resp,cuesampr,buttonpress,data.fsample);
            
            
            flan = 'LI';
            resp = 'right';
            % buttonpress to cue
            cuesampl = trig_sample.sample(strcmp(trig_sample.type, flan) ) ;
            [respLIcorr,flanLIcorr] = resp2trial(resp,cuesampl,buttonpress,data.fsample);
            
            resp = 'left';
            [respLIcomm,flanLIcomm] = resp2trial(resp,cuesampl,buttonpress,data.fsample);
            
            
            flan = 'LC';
            resp = 'left';
            % buttonpress to cue
            cuesamp = trig_sample.sample(strcmp(trig_sample.type, flan) ) ;
            [respLCcorr,flanLCcorr] = resp2trial(resp,cuesamp,buttonpress,data.fsample);
            
            flan = 'RC';
            resp = 'right';
            % buttonpress to cue
            cuesamp = trig_sample.sample(strcmp(trig_sample.type, flan) ) ;
            [respRCcorr,flanRCcorr] = resp2trial(resp,cuesamp,buttonpress,data.fsample);
            
            flan = 'LC';
            resp = 'right';
            % buttonpress to cue
            cuesamp = trig_sample.sample(strcmp(trig_sample.type, flan) ) ;
            [respLCcomm,flanLCcomm] = resp2trial(resp,cuesamp,buttonpress,data.fsample);
            
            flan = 'RC';
            resp = 'left';
            % buttonpress to cue
            cuesamp = trig_sample.sample(strcmp(trig_sample.type, flan) ) ;
            [respRCcomm,flanRCcomm] = resp2trial(resp,cuesamp,buttonpress,data.fsample);
            %% Co-register MRI
            
            
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
            
            
            %% Theta power
            filt_order = []; % default
            % data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,50,filt_order,'but');
            data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[4 8],filt_order,'but');
            dataf = data;
            dataf.trial{1} = data_filt;
            
            timewflan =[0.25 0.65];
            timewbase =[-0.45 -0.05];
            datai = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],dataf,time,timewflan,1);
            Ci = cov(cell2mat(datai.trial)');
            
            dataib = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],dataf,time,timewbase,1);
            Cib = cov(cell2mat(dataib.trial)');
            
            dataC = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],dataf,time,timewflan,1);
            Cc = cov(cell2mat(dataC.trial)');
            
            dataCb = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],dataf,time,timewbase,1);
            Ccb = cov(cell2mat(dataCb.trial)');
            
            
            timewresp = [0, 0.5];
            dataA = define_trials([respRIcomm;respLIcomm],dataf,time,timewresp,1);
            Cerr = cov(cell2mat(dataA.trial)');
            
            dataC = define_trials([respRIcorr;respLIcorr],dataf,time,timewresp,1);
            Ccorr = cov(cell2mat(dataC.trial)');
            
            P = cell(size(L));
            Pc = cell(size(L));
            Pi = cell(size(L));
            Perr = cell(size(L));
            
            for l = 1:length(W)
                
                w = W{l};
                
                P{l} = ((w'*Ci*w) - (w'*Cc*w) ) / (2*w'*noiseCtheta*w);
                Pc{l} = ((w'*Cc*w) - (w'*Ccb*w) ) / (2*w'*noiseCtheta*w);
                Pi{l} = ((w'*Ci*w) - (w'*Cib*w) ) / (2*w'*noiseCtheta*w);
                Perr{l} = ((w'*Cerr*w) - (w'*Ccorr*w) ) / (2*w'*noiseCtheta*w);
                %     P{l} = (w'*Ca*w)  / (w'*noiseC*w);
                %     Pc{ii} = (w'*Cc*w)  / (w'*noiseC*w);
                
                
            end
            
            Theta_Flan_Inc_Cong  = cell2mat(P');
            Theta_Flan_Inc  = cell2mat(Pi');
            Theta_Flan_Cong  = cell2mat(Pc');
            Theta_Resp_Err_Corr  = cell2mat(Perr');
            
            theta_time = timewflan;
            
            
            % Alpha
            filt_order = []; % default
            % data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,50,filt_order,'but');
            data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[8 13],filt_order,'but');
            dataf = data;
            dataf.trial{1} = data_filt;
            
            timewflan =[0.25 0.65];
            timewbase =[-0.45 -0.05];
            datai = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],dataf,time,timewflan,1);
            Ci = cov(cell2mat(datai.trial)');
            
            dataib = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],dataf,time,timewbase,1);
            Cib = cov(cell2mat(dataib.trial)');
            
            dataC = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],dataf,time,timewflan,1);
            Cc = cov(cell2mat(dataC.trial)');
            
            dataCb = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],dataf,time,timewbase,1);
            Ccb = cov(cell2mat(dataCb.trial)');
            
            
            P = cell(size(L));
            Pc = cell(size(L));
            Pi = cell(size(L));
            
            for l = 1:length(W)
                
                w = W{l};
                
                P{l} = ((w'*Ci*w) - (w'*Cc*w) ) / (2*w'*noiseCalpha*w);
                Pc{l} = ((w'*Cc*w) - (w'*Ccb*w) ) / (2*w'*noiseCalpha*w);
                Pi{l} = ((w'*Ci*w) - (w'*Cib*w) ) / (2*w'*noiseCalpha*w);
                %     P{l} = (w'*Ca*w)  / (w'*noiseC*w);
                %     Pc{ii} = (w'*Cc*w)  / (w'*noiseC*w);
                
                
            end
            
            Alpha_Flan_Inc_Cong  = cell2mat(P');
            Alpha_Flan_Cong  = cell2mat(Pc');
            Alpha_Flan_Inc  = cell2mat(Pi');
            
            alpha_time = timewresp;
            
            
            save( sprintf('%s%s/Congruency_multiSpheres_%dmm_regmu%d.mat',...
                processingfolder,filename(1:end-3),gridres,mu) ,...
                'Theta_Flan_Inc_Cong','Theta_Flan_Cong','Theta_Flan_Inc','theta_time','Theta_Resp_Err_Corr',...
                'Alpha_Flan_Inc_Cong','Alpha_Flan_Cong','Alpha_Flan_Inc','alpha_time');
            
            
            %% ERN localization
            filt_order = []; % default
            data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,50,filt_order,'but');
            dataf = data;
            dataf.trial{1} = data_filt;
            
            
            timewresp =[-0.5 0.8];
            dataA = define_trials([respRCcomm;respLCcomm;respRIcomm;respLIcomm],dataf,time,timewresp,1);
            
            dataC = define_trials([respRCcorr;respLCcorr],dataf,time,timewresp,1);
            dataI = define_trials([respRIcorr;respLIcorr],dataf,time,timewresp,1);
            
            cfg = [];
            cfg.keeptrials = 'no';
            erpCcorr = ft_timelockanalysis([], dataC);
            erpIcorr = ft_timelockanalysis([], dataI);
            erpcomm = ft_timelockanalysis([], dataA);
            
            erpcommav=  erpcomm.avg;
            erpCcorrav =  erpCcorr.avg;
            erpIcorrav =  erpIcorr.avg;
            
            Pcomm = cell(size(L));
            PCcorr = cell(size(L));
            PIcorr = cell(size(L));
            
            for l = 1:length(W)
                w = Wdc{l};% needs depth corrected weights
                
                Pcomm{l} = (w'*erpcommav);
                PCcorr{l} = (w'*erpCcorrav);
                PIcorr{l} = (w'*erpIcorrav);
            end
            
            Pcomm  = cell2mat(Pcomm');
            PCcorr  = cell2mat(PCcorr');
            PIcorr  = cell2mat(PIcorr');
            
            resptime = erpcomm.time;
            
            
            %% Flanker ERF localization
            
            timewresp =[-0.2 0.8];
            
            dataC = define_trials([flanRCcorr;flanLCcorr;flanRCcomm;flanLCcomm],dataf,time,timewresp,1);
            dataI = define_trials([flanRIcorr;flanLIcorr;flanRIcomm;flanLIcomm],dataf,time,timewresp,1);
            
            cfg = [];
            cfg.keeptrials = 'no';
            erpCcorr = ft_timelockanalysis([], dataC);
            erpIcorr = ft_timelockanalysis([], dataI);
            
            erpCcorrav =  erpCcorr.avg;
            erpIcorrav =  erpIcorr.avg;
            
            PCflan = cell(size(L));
            PIflan = cell(size(L));
            
            for l = 1:length(W)
                w = Wdc{l};% needs depth corrected weights
                
                PCflan{l} = (w'*erpCcorrav);
                PIflan{l} = (w'*erpIcorrav);
            end
            
            PCflan  = cell2mat(PCflan');
            PIflan  = cell2mat(PIflan');
            
            flantime = erpCcorr.time;
            save( sprintf('%s%s/ERN_multiSpheres_%dmm_regmu%d.mat',...
                processingfolder,filename(1:end-3),gridres,mu) ,'Pcomm','PIcorr','PCcorr','PCflan','PIflan','flantime','resptime' );
            
        end
    end
end
%%

