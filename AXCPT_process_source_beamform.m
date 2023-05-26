%function AXCPT_process_source_beamform(filenames,datapath,processingfolder,bvfiles,fids_name)

addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

sublist = {'24531';'24563';'24590';'24580';'24626';'24581';'24482';...
    '24640';'24592';'24667';'24678' };

for ss = 1:length(sublist)
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
        if  contains(d(ii).name,'meg.ds') && contains(d(ii).name,'task-axcpt')
            n = n+1;
            filenames{n} = d(ii).name;
        end
       
    end
   
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
    
    d = dir('/data/EDB/MEG_AXCPT_Flanker/data/emptyroom/');
    emptyroom = []; jj = 0;
    for ii = 3:length(d)
        if contains(d(ii).name, TaskDate(1:8))
            jj = jj + 1;
            emptyroom{jj} = ['/data/EDB/MEG_AXCPT_Flanker/data/emptyroom/',d(ii).name];
        end
    end
    
    
    highpass = 1;
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
    % noiseCtheta = zeros(length(channels),length(channels),length(emptyroom));
    for ii = 1:length(emptyroom)
        
        
        cfg = [];
        cfg.dataset = emptyroom{ii};
        cfg.continuous = 'yes';
        cfg.channel = channels;
        cfg.demean = 'yes';
        cfg.detrend = 'no';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [lowpass, highpass]; % With notch filter 60Hz
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [58 62; 118 122; 178 182]; % With notch filter 60Hz
        
        data_empty = ft_preprocessing(cfg);
        emptyC = zeros(length(channels),length(channels),length(data_empty.trial));
        for t = 1:length(data_empty.trial)
            emptyC = cov(data_empty.trial{t}');
        end
        noiseC(:,:,ii) = mean(emptyC,3);
        
        
        
        %     cfg.bpfreq = [4, 8]; % With notch filter 60Hz
        %     data_empty = ft_preprocessing(cfg);
        %     emptyC = zeros(length(channels),length(channels),length(data_empty.trial));
        %     for t = 1:length(data_empty.trial)
        %         emptyC = cov(data_empty.trial{t}');
        %     end
        %     noiseCtheta(:,:,ii) = mean(emptyC,3);
        
    end
    noiseC = mean(noiseC,3);
    % noiseCtheta = mean(noiseCtheta,3);
    
    
    %% Standard pre-processing
    
    cd(datapath)
    
    if ~exist(mri_name,'file')
        mri_name = [mri_name,'.gz'];
    end
    fids_file =  [datapath(1:end-4),'anat/',fids_name];
    if ~exist(fids_file,'file')
        fids_file =  [datapath(1:end-4),'anat/sub-',sub,'_fiducials_axcpt.tag'];
    end
    mri = fids2ctf(mri_name,fids_file,0);

    
    for ii = 1:length(filenames)
        filename = filenames{ii};
        sub = filename(5:9);
        
        %     if ~exist([processingfolder,'/',filename(1:end-3),'/ICA_artifacts.mat'], 'file')
        
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
        
        
        %% Find samples
        cue = 'A';
        % buttonpress to cue
        cuesampA = trig_sample.sample(strcmp(trig_sample.type, cue) ) ;
        
        [buttonsampA,rtA] = match_responses(cuesampA, buttonpress, 'left', data.fsample);
        cuesampA(buttonsampA==0) = [];
        
        cue = 'B';
        % buttonpress to cue
        cuesampB = trig_sample.sample(strcmp(trig_sample.type, cue) ) ;
        [buttonsampB,rtB] = match_responses(cuesampB, buttonpress, 'left', data.fsample);
        cuesampB(buttonsampB==0) = [];
        
        probe = 'AY';
        %buttonpress to probe
        probesampAY = trig_sample.sample(strcmp(trig_sample.type, probe)) ;
        
        buttonsampAYcomm = match_responses(probesampAY, buttonpress, 'right', data.fsample);
        [buttonsampAYcorr, rtAY] = match_responses(probesampAY, buttonpress, 'left', data.fsample);
        probesampAY( buttonsampAYcomm == 0 & buttonsampAYcorr == 0 ) = [];
        
        
        probe = 'AX';
        %buttonpress to probe
        probesampAX = trig_sample.sample(strcmp(trig_sample.type, probe)) ;
        
        buttonsampAXcorr = match_responses(probesampAX, buttonpress, 'right', data.fsample);
        buttonsampAXcomm = match_responses(probesampAX, buttonpress, 'left', data.fsample);
        probesampAX( buttonsampAXcomm == 0 & buttonsampAXcorr == 0 ) = [];
        
        probe = 'BX';
        %buttonpress to probe
        probesampBX = trig_sample.sample(strcmp(trig_sample.type, probe)) ;
        
        buttonsampBXcomm = match_responses(probesampBX, buttonpress, 'right', data.fsample);
        buttonsampBXcorr = match_responses(probesampBX, buttonpress, 'left', data.fsample);
        probesampBX( buttonsampBXcomm == 0 & buttonsampBXcorr == 0 ) = [];
        
        probe = 'BY';
        %buttonpress to probe
        probesampBY = trig_sample.sample(strcmp(trig_sample.type, probe)) ;
        
        buttonsampBYcomm = match_responses(probesampBY, buttonpress, 'right', data.fsample);
        buttonsampBYcorr = match_responses(probesampBY, buttonpress, 'left', data.fsample);
        probesampBY( buttonsampBYcomm == 0 & buttonsampBYcorr == 0 ) = [];
        
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
        
        
        %% ERF
        freq = 50;
        filt_order = []; % default
        
        dataf = data;
        data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,freq,filt_order,'but');
        dataf.trial{1} = data_filt;
        clear data_filt
        
        twind = [-0.5 1];
        [dataprobeAY,ttdel] = define_trials(probesampAY ,dataf,time,twind,1);
        [dataprobeAX,ttdel] = define_trials(probesampAX,dataf,time,twind,1);
        [dataprobeBY,ttdel] = define_trials(probesampBY ,dataf,time,twind,1);
        [dataprobeBX,ttdel] = define_trials(probesampBX,dataf,time,twind,1);
        [datacueA,ttdel] = define_trials(cuesampA,dataf,time,twind,1);
        [datacueB,ttdel] = define_trials(cuesampB,dataf,time,twind,1);
        
        ntrials = [];
        ntrials.A = length(datacueA.trial);
        ntrials.B = length(datacueB.trial);
        ntrials.AX = length(dataprobeAX.trial);
        ntrials.AY = length(dataprobeAY.trial);
        ntrials.BX = length(dataprobeBX.trial);
        ntrials.BY = length(dataprobeBY.trial);
        
        %     [datapressAY,ttdel] = define_trials([buttonsampAYcorr ],dataf,time,twind,1);
        %     [datapressAX,ttdel] = define_trials([buttonsampAXcorr ],dataf,time,twind,1);
        %     [datapressBY,ttdel] = define_trials( [buttonsampBYcorr ],dataf,time,twind,1);
        %     [datapressBX,ttdel] = define_trials([buttonsampBXcorr ],dataf,time,twind,1);
        %     [datapressA,ttdel] = define_trials(buttonsampA,dataf,time,twind,1);
        %     [datapressB,ttdel] = define_trials(buttonsampB,dataf,time,twind,1);
        %
        %     [datapresserr,ttdel] = define_trials([buttonsampAYcomm; buttonsampAXcomm;buttonsampBXcomm;buttonsampBYcomm],dataf,time,twind,1);
        
        cfg = [];
        cfg.demean        = 'yes';
        cfg.baselinewindow = [-inf, -0];
        
        if strcmp(cfg.demean, 'yes')
            datacueA = ft_preprocessing(cfg,datacueA);
            datacueB = ft_preprocessing(cfg,datacueB);
            dataprobeAX = ft_preprocessing(cfg,dataprobeAX);
            dataprobeAY = ft_preprocessing(cfg,dataprobeAY);
            dataprobeBX = ft_preprocessing(cfg,dataprobeBX);
            dataprobeBY = ft_preprocessing(cfg,dataprobeBY);
            
            %         datapressA = ft_preprocessing(cfg,datapressA);
            %         datapressB = ft_preprocessing(cfg,datapressB);
            %         datapressAX = ft_preprocessing(cfg,datapressAX);
            %         datapressAY = ft_preprocessing(cfg,datapressAY);
            %         datapressBX = ft_preprocessing(cfg,datapressBX);
            %         datapressBY = ft_preprocessing(cfg,datapressBY);
            
            %         datapresserr = ft_preprocessing(cfg,datapresserr);
            
        end
        
        datacueA = ft_timelockanalysis([], datacueA);
        datacueB = ft_timelockanalysis([],datacueB );
        dataprobeAX = ft_timelockanalysis([], dataprobeAX);
        dataprobeAY = ft_timelockanalysis([],dataprobeAY );
        dataprobeBX = ft_timelockanalysis([], dataprobeBX);
        dataprobeBY = ft_timelockanalysis([],dataprobeBY );
        
        
        %     datapressA = ft_timelockanalysis([], datapressA);
        %     datapressB = ft_timelockanalysis([],datapressB );
        %     datapressAX = ft_timelockanalysis([], datapressAX);
        %     datapressAY = ft_timelockanalysis([],datapressAY );
        %     datapressBX = ft_timelockanalysis([], datapressBX);
        %     datapressBY = ft_timelockanalysis([],datapressBY );
        
        
        %     datapresserr = ft_timelockanalysis([],datapresserr );
        
        erpcueA=  datacueA.avg;
        erpcueB =  datacueB.avg;
        erpprobeAX =  dataprobeAX.avg;
        erpprobeAY =  dataprobeAY.avg;
        erpprobeBX =  dataprobeBX.avg;
        erpprobeBY =  dataprobeBY.avg;
        
        %     erppressA=  datapressA.avg;
        %     erppressB =  datapressB.avg;
        %     erppressAX =  datapressAX.avg;
        %     erppressAY =  datapressAY.avg;
        %     erppressBX =  datapressBX.avg;
        %     erppressBY =  datapressBY.avg;
        
        %     erppresserr =  datapresserr.avg;
        
        
        PerpAX  = cell(size(L));
        PerpAY  = cell(size(L));
        PerpBX  = cell(size(L));
        PerpBY  = cell(size(L));
        PerpA  = cell(size(L));
        PerpB  =cell(size(L));
        
        %     BerpAX  = cell(size(L));
        %     BerpAY  = cell(size(L));
        %     BerpBX  = cell(size(L));
        %     BerpBY  = cell(size(L));
        %     BerpA  = cell(size(L));
        %     BerpB  =cell(size(L));
        
        %     Berperr  = cell(size(L));
        
        for l = 1:length(W)
            w = Wdc{l};% needs depth corrected weights
            
            PerpAX{l} = (w'*erpprobeAX);
            PerpAY{l} = (w'*erpprobeAY);
            PerpBX{l} = (w'*erpprobeBX);
            PerpBY{l} = (w'*erpprobeBY);
            PerpA{l} = (w'*erpcueA);
            PerpB{l} = (w'*erpcueB);
            
            %         BerpAX{l} = (w'*erppressAX);
            %         BerpAY{l} = (w'*erppressAY);
            %         BerpBX{l} = (w'*erppressBX);
            %         BerpBY{l} = (w'*erppressBY);
            %         BerpA{l} = (w'*erppressA);
            %         BerpB{l} = (w'*erppressB);
            
            %         Berperr{l} = (w'*erppresserr);
            
            
        end
        
        PerpAX  = cell2mat(PerpAX');
        PerpAY  = cell2mat(PerpAY');
        PerpBX = cell2mat(PerpBX');
        PerpBY = cell2mat(PerpBY');
        PerpA  = cell2mat(PerpA');
        PerpB  = cell2mat(PerpB');
        
        %     BerpAX  = cell2mat(BerpAX');
        %     BerpAY  = cell2mat(BerpAY');
        %     BerpBX = cell2mat(BerpBX');
        %     BerpBY = cell2mat(BerpBY');
        %     BerpA  = cell2mat(BerpA');
        %     BerpB  = cell2mat(BerpB');
        
        %     Berperr  = cell2mat(Berperr');
        
        erptime = datacueA.time;
        
        %     save( sprintf('%s%s/ERF_multiSpheres_%dmm_regmu%d.mat',...
        %         processingfolder,filename(1:end-3),gridres,mu) ,...
        %         'PerpAX','PerpAY','PerpBX','PerpBY','PerpA','PerpB','erptime',...
        %         'BerpAX','BerpAY','BerpBX','BerpBY','BerpA','BerpB','ntrials');
        
        save( sprintf('%s%s/ERF_multiSpheres_%dmm_regmu%d.mat',...
            processingfolder,filename(1:end-3),gridres,mu) ,...
            'PerpAX','PerpAY','PerpBX','PerpBY','PerpA','PerpB','erptime','ntrials');
        %% Theta power
        %
        %     freq = [4 8];
        %     filt_order = []; % default
        %
        %     dataf = data;
        %     data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,freq,filt_order,'but');
        %     dataf.trial{1} = data_filt;
        %     clear data_filt
        %
        %     twind = [0.1 0.5];
        %     [dataprobeAY,~] = define_trials(probesampAY ,dataf,time,twind,1);
        %     [dataprobeAX,~] = define_trials(probesampAX,dataf,time,twind,1);
        %     [dataprobeBY,~] = define_trials(probesampBY ,dataf,time,twind,1);
        %     [dataprobeBX,~] = define_trials(probesampBX,dataf,time,twind,1);
        %     [datacueA,~] = define_trials(cuesampA,dataf,time,twind,1);
        %     [datacueB,~] = define_trials(cuesampB,dataf,time,twind,1);
        %
        %     CA = cov(cell2mat(datacueA.trial)');
        %     CB = cov(cell2mat(datacueB.trial)');
        %     CAX = cov(cell2mat(dataprobeAX.trial)');
        %     CAY = cov(cell2mat(dataprobeAY.trial)');
        %     CBX = cov(cell2mat(dataprobeBX.trial)');
        %     CBY = cov(cell2mat(dataprobeBY.trial)');
        %
        %
        %     PprobeAX = cell(size(L));
        %     PprobeAY = cell(size(L));
        %     PprobeBX = cell(size(L));
        %     PprobeBY = cell(size(L));
        %     PcueA = cell(size(L));
        %     PcueB = cell(size(L));
        %     for l = 1:length(L)
        %         w = W{l};
        %         PprobeAY{l} = (w'*CAY*w)  / (w'*noiseCtheta*w);
        %         PprobeAX{l} = (w'*CAX*w)  / (w'*noiseCtheta*w);
        %         PprobeBY{l} = (w'*CBY*w)  / (w'*noiseCtheta*w);
        %         PprobeBX{l} = (w'*CBX*w)  / (w'*noiseCtheta*w);
        %         PcueA{l} = (w'*CA*w) / (w'*noiseCtheta*w);
        %         PcueB{l} = (w'*CB*w) / (w'*noiseCtheta*w);
        %     end
        %
        %     PprobeAX  = cell2mat(PprobeAX)';
        %     PprobeAY  = cell2mat(PprobeAY)';
        %     PprobeBX  = cell2mat(PprobeBX)';
        %     PprobeBY  = cell2mat(PprobeBY)';
        %     PcueA  = cell2mat(PcueA)';
        %     PcueB  = cell2mat(PcueB)';
        %
        %     save( sprintf('%s%s/Theta_multiSpheres_%dmm_regmu%d.mat',...
        %         processingfolder,filename(1:end-3),gridres,mu) ,...
        %         'PprobeAX','PprobeAY','PprobeBX','PprobeBY','PcueA','PcueB');
        %     end
    end
end