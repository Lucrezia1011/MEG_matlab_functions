function [data,BadSamplesAll] = preproc_bids_2024(data_name,data_path,highpass,lowpass,icaopt,plotopt)
% Pre-processing 
% [data,BadSamplesAll] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt)
% data_name     = dataset name (.ds)
% highpass      = highpass frequency in Hz
% lowpass       = lowpass frequency in Hz
% icaopt        = 1 to use ICA to remove eyeblinks and heartbeat
% plotopt       = 1 to plot processed MEG data
%
% Warning: hard-coded derivates data_path
addpath('/home/liuzzil2/fieldtrip-20190812/fieldtrip_private')

% Data header info
hdr = ft_read_header(data_name);
% Get Bad channel names
fid = fopen([data_name,'/BadChannels']);
BadChannels = textscan(fid,'%s');
fclose(fid);

% get MEG channel names
channels = hdr.label(strncmp(hdr.label,'M',1));
% Delete Bad channels
chanInd = zeros(size(channels));
for iiC = 1:length(BadChannels{1})
    chanInd = chanInd | strcmp(channels,BadChannels{1}{iiC});
end
channels(find(chanInd)) = [];

%% Finds bad data segments
 
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = channels;
cfg.demean = 'yes';
cfg.detrend = 'no';
cfg.bpfilter = 'no';
cfg.bsfilter = 'yes';
cfg.bsfreq = [58 62; 118 122; 178 182]; % With notch filter 60Hz

data = ft_preprocessing(cfg);

f = data.fsample;
% Eliminate bad trials(extract from ft_read_event)
[condNumbers,condLabels] = read_ctf_cls([data_name,'/ClassFile.cls']);

if any(strcmp(condLabels,'BAD'))
    condBad = condNumbers{strcmp(condLabels,'BAD')};
else
    condBad = [];
end
if any(strcmp(condLabels,'BAD_HeadMotion'))
    condBadHead =  condNumbers{strcmp(condLabels,'BAD_HeadMotion')};
else
    condBadHead = [];
end
% Find Trials marked as BAD
BadNumbers = unique([condBad, condBadHead]);

if ~isempty(BadNumbers) 
    BadTrials = zeros(data.hdr.nSamples,length(BadNumbers));
    for jN=1:length(BadNumbers)        
        BadTrials(:,jN) = (BadNumbers(jN)-1)*data.hdr.nSamples + (1:data.hdr.nSamples) ;       
    end
    BadTrials = BadTrials(:);
else
    BadTrials = [];
end

% Find Bad data segments 
bsegFile =  [data_name,'/bad.segments'];
if exist(bsegFile,'file')

    fid = fopen(bsegFile);
    bsegs = (textscan(fid,'%s'));
    fclose(fid);
    bsegs = bsegs{1};
      % delete fit error
    ferr =  find(strcmp(bsegs,'FitError') );
    fdel = zeros(length(ferr),4);
    for fr = 1:length(ferr)
        fdel(fr,:) = ferr(fr) + [-3:0];
    end
    bsegs(fdel(:)) = [];

    % find headmotion labels
    ferr =  find(strcmp(bsegs,'HeadMotion') );
    BadMotion = cell(size(ferr,1),1);
    if ~isempty(ferr)   
        for jN=1:length(BadMotion)    
            bsegM = str2double(bsegs(ferr(jN) + [-3:-1]));
            BadMotion{jN} = (bsegM(1)-1)*data.hdr.nSamples +...
                (round(bsegM(2)*f):round(bsegM(3)*f)) +1 ;       
        end  
    end
    BadMotionall = cell2mat(BadMotion');
    for jN = 1:length(BadMotion)
        bmw = floor(median(BadMotion{jN}))+(-10*f/2:10*f/2); % look into 10s windows
        bmb = intersect(bmw,BadMotionall); % section labelled with bad motion
        if nnz(bmb) > (10*f/2) % if at least 5s into a continus 10s window is labelled as bad motion
            BadMotion{jN} = union(bmw,BadMotion{jN});
        else
            BadMotion{jN} = [];
        end
    end
    BadMotion = cell2mat(BadMotion');
    BadMotion(BadMotion > data.sampleinfo(2)) = [];
    % delete head motion from bad segments
    ferr =  find(strcmp(bsegs,'HeadMotion') );
    fdel = zeros(length(ferr),4);
    for fr = 1:length(ferr)
        fdel(fr,:) = ferr(fr) + [-3:0];
        
    end
    bsegs(fdel(:)) = [];
   

    bsegs = reshape(str2double(bsegs),3,size(bsegs,1)/3)';
    BadSegs = cell(size(bsegs,1),1);
    if ~isempty(bsegs)   
        for jN=1:length(BadSegs)        
            BadSegs{jN} = (bsegs(jN,1)-1)*data.hdr.nSamples +...
                (round(bsegs(jN,2)*f):round(bsegs(jN,3)*f)) +1 ;       
        end  
    end
    BadSegs = cell2mat(BadSegs');
else
    BadSegs = [];
    BadMotion =[];
end
BadEdges = [1:f, data.sampleinfo(2)+(-f:0)]'; % delete edges after notch filter
% Combine Bad Trials and Bad Segments
BadSamplesAll = unique([BadTrials; BadSegs';BadMotion';BadEdges]);
% Find Bad tail of dataset
indLast = find(diff(BadSamplesAll)~=1);
BadSamplesCell = cell(nnz(indLast)+1,1);

if ~isempty(BadSamplesAll)
    for jN = 1:nnz(indLast)+1
        if jN == 1 && isempty(indLast)
            BadSamplesCell{jN} = BadSamplesAll;
        elseif jN == 1 && ~isempty(indLast)
             BadSamplesCell{jN} = BadSamplesAll(1):BadSamplesAll(indLast(1));
        elseif jN == length(indLast)+1  && ~isempty(indLast)
            BadSamplesCell{jN} = BadSamplesAll(indLast(jN-1)+1):BadSamplesAll(end);
        else
            BadSamplesCell{jN} = BadSamplesAll(indLast(jN-1)+1):BadSamplesAll(indLast(jN));
        end
        
    end
end


if ~isempty(BadSamplesAll)
    if BadSamplesAll(end) == data.sampleinfo(2)
        BadSamplesLast = BadSamplesCell{end};
        BadSamplesCell(end) = [];
        BadSamples = cell2mat(BadSamplesCell');
    
        % Eliminate Bad Segments at end of dataset including aborted data
        % collection
        data.time{1}(BadSamplesLast) = [];
        data.trial{1}(:,BadSamplesLast) = [];
        data.sampleinfo = [1 length(data.time{1})];
%         indLast = find(diff(BadSamples)~=1);
    else
        BadSamples = BadSamplesAll;
    end
   
    
else
    BadSamples = BadSamplesAll;
end

% df = diff(data.trial{1},1,2);
% figure; histogram(df(:))
% % distribution of consecutive time points has stdev = 73.17 fT

%% Filter and eliminate bad segments


% Find large SQUID jumps
[sampleJump,sensJump] = find(abs(diff(data.trial{1}'))>20e-12); 


if ~isempty(sensJump)

    BadSQUID = cell(0);
    iiN = 0;
    for iib = 1:length(BadSamplesCell)
    
        iia = ismember(sampleJump,BadSamplesCell{iib});
        if nnz(iia) > 0
            iiN = iiN+1;
            BadSQUID{iiN} = BadSamplesCell{iib};
        end
        
    end
    if iiN == 0 
        warning('Sensors with jumps: ')
        disp(channels(unique(sensJump)))
        warning('Time of jumps: ')
        disp(channels(unique(sampleJump/f)))
        error('Found large data jumps unmarked in raw data')
    end

    dataBadSQUID = data.trial{1};

    % correct jumps
    for iS = 1:length(BadSQUID)
      
        s2 = BadSQUID{iS}(1)-1;
        % mean correct
        if s2 == 0
            m1 = 0;
        else
            m1 = mean(dataBadSQUID(:,(-3*f:0)+s2),2);
        end
%         dataBadSQUID(:,s1:s2) = dataBadSQUID(:,s1:s2)-m1;

        s1 = BadSQUID{iS}(end)+1;
        s2 = data.sampleinfo(2);
        % mean correct
        m2 = mean(dataBadSQUID(:,(0:3*f)+s1),2);
        dataBadSQUID(:,s1:s2) = dataBadSQUID(:,s1:s2)-m2+m1;

%         BadSamples(BadSamples == BadSQUID{iS}) = [];

    end
   
    dataBadSQUID(:,cell2mat(BadSQUID)) =[];
    dataBadSQUID = dataBadSQUID - mean(dataBadSQUID,2);
    filt_order = 4; % default
    if highpass <=1 && lowpass <= 30
        filt_order = 3;
    end
    dataBadSQUID = ft_preproc_bandpassfilter(dataBadSQUID, f, [highpass lowpass], filt_order, [], [], []);

    
    % pad data channels with deleted jumps and substite into data struct
    indN = true(1,data.sampleinfo(2));
    for iS = 1:length(BadSQUID)
        indN(BadSQUID{iS}) = false;
    end
  
    data.trial{1}(:,indN) = dataBadSQUID;
    data.trial{1}(:,~indN) = 0;
    
    warning('Found %d SQUID jumps.\n',iiN)
   
    clear dataBadSQUID 
else


    % Filter data
    data.trial{1} = data.trial{1} - mean(data.trial{1},2);
    % Are we introducing filtering artefacts by applying highpass filter on
    % data with discontinuities (i.e. after deleting bad segments)? Yes
    % Apply highpass filter after correcting jumps, but before introducing
    % discontinuities from bad data segments
    filt_order = 4; % default
    if highpass<=1 && lowpass <= 30
        filt_order = 3;
    end
    data.trial{1} = ft_preproc_bandpassfilter(data.trial{1}, f, [highpass lowpass], filt_order, [], [], []);

end
% Outlier identification?
% stdev = std(data.trial{1},[],2);
% stdevt = abs(data.trial{1}) > stdev*10;
% [ch,ss] = ind2sub(size(stdevt), find(stdevt));
figure(1); clf; 
if plotopt == 1
    set(gcf,'Position',[571   231   577   705])
    subplot(311)
    plot(data.trial{1}')
    title('MEG data with bad segments')
end

% Eliminate Bad Segments
BadSamples(BadSamples<1) = [];
BadSamplesAll(BadSamplesAll<1) = [];
data.time{1}(BadSamples) = [];
data.trial{1}(:,BadSamples) = [];
data.sampleinfo = [1 length(data.time{1})];

if plotopt == 1
     subplot(312)
    plot(data.trial{1}')
    title('MEG data after bad segments deletion')
    drawnow
end

%% ICA
% if ~isempty(sampleJump) 
%     keyboard %dbcont
% end
if icaopt == 1
   
    data_path = [data_path,'/',data_name(1:end-3)];
    if ~exist(data_path,'dir')
        mkdir(data_path)
    end
    
    rmpath(('~/fieldtrip-20190812/fieldtrip_private'))
    if exist([data_path,'/ICA_artifacts.mat'],'file')
        load([data_path,'/ICA_artifacts.mat']);
    else
        if any(strcmp(data.hdr.label,'UADC009'))
            cfg = [];
            cfg.dataset = data_name;
            cfg.continuous = 'yes';
            cfg.channel = {'UADC009';'UADC010';'UADC013'};
            eyed = ft_preprocessing(cfg);
            eye = eyed.trial{1};
        end
        dataica = data;
        if highpass < 1
            dataica.trial{1} = ft_preproc_highpassfilter(data.trial{1}, f, 1, filt_order, [], [], []);     
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
        elseif icomp < 20
            icomp = 20;
        end
        cfg =[];
        
        cfg.method = 'fastica';
        cfg.fastica.numOfIC = icomp;
        cfg.fastica.maxNumIterations = 100;
        comp = ft_componentanalysis(cfg, dataica);
        icomp = length(comp.label);
        
        figure(2)
        cfg           = [];
        cfg.component = [1:icomp];       % specify the component(s) that should be plotted
        cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)
        set(gcf,'position',[ 61    76   964   847])
        
        cfg          = [];
        cfg.channel  = [1:5]; % components to be plotted
        cfg.viewmode = 'component';
        cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
        cfg.continuous   = 'yes' ;
        cfg.blocksize    = 10;
        ft_databrowser(cfg, comp)
        set(gcf,'position',[ 571         198        1030         742])

        if exist('eye','var')
            figure(4); % delete bad data segments
            eye(:,BadSamplesAll) = [];
            plot(abs(corr(eye',comp.trial{1}'))','*')
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
        save([data_path,'/ICA_artifacts'],'comps')
        
        close([2,3,4])
        
        
    end
    
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data          = ft_rejectcomponent(cfg, comps,data);
    if plotopt == 1
        figure(1); subplot(313);plot(data.trial{1}')
        title('MEG data after ICA'); drawnow

    else
        close all
    end
end

% Delete Bad channels from gradiometers structure
chanInd = zeros(size(data.grad.label));
for iiC = 1:length(BadChannels{1})
    chanInd = chanInd | strcmp(data.grad.label,BadChannels{1}{iiC});
end
data.grad.label(find(chanInd)) = [];
data.grad.chanori(find(chanInd),:) = [];
data.grad.chanpos(find(chanInd),:) = [];
data.grad.tra(find(chanInd),:) = [];
