function [samples_T,trig_sample,buttonpress] = matchTriggers( filename, BadSamples)
% Lucrezia Liuzzi, last updated 2024/11/01
%
% Matches accurate time of projector display to data triggers
%
% [samples_T,trig_sample,buttonpress] = matchTriggers(data_name, BadSamples)
% filename       = name of task-mmi3 dataset (.ds)
% BadSamples    = time samples to exclude, output of 'preproc_bids'
% samples_T     = matching samples after deleting bad samples
% trig_sample   = output structure with all triggers with original samples
% buttonpress   = output button presses



%% Select data

% Outdated way to read light pixel channel
% pix = ft_read_data(data_name,'chanindx','UADC016');

cfg = [];
cfg.dataset = filename;
cfg.continuous = 'no';
cfg.channel = 'SCLK01'; % Time channel!
time = ft_preprocessing(cfg);


% Samples transform from continuous recording to pre-processed (no Bad segs)
samples_cut = 1:(time.sampleinfo(end)-nnz(BadSamples));
indN = true(1,time.sampleinfo(end));
indN(BadSamples)  =false;
samples_T = zeros(time.sampleinfo(end),1);
samples_T(indN) = samples_cut;

event = ft_read_event(filename);
if size(event,1) == 1
    event = event';
end
% read Trigger channel in MEG data

trig_sample = struct;
trig_sample.sample = zeros(size(event,1),1);
trig_sample.value = zeros(size(event,1),1);
trig_sample.type = cell(size(event,1),1);
ii = 0;
if any(contains({event.type},'UPPT001'))
    for tt = 1:length(event)
        if strcmp(event(tt).type, 'UPPT001' )
            ii = ii +1;
            trig_sample.sample(ii,1) = event(tt).sample ; % Average 24 sample delay
            trig_sample.value(ii,1) = event(tt).value;
            trig_sample.type{ii,1} = event(tt+1).type;
        end
    end
    trig_sample.sample = trig_sample.sample(1:ii,:);
    trig_sample.value = trig_sample.value(1:ii,:);
    trig_sample.type = trig_sample.type(1:ii,:);
else
    event(contains({event.type},'bad')) = []; % with ft2024 dimensions need to be flipped
    trig_sample.sample = [event.sample]';
    trig_sample.value = [event.value]';
    trig_sample.type = {event.type}';
end

if ~any(contains(trig_sample.type,'resp_')) % already adjusted times from bids conversion
    
    % read LIGHT marker (projector display time)
    if any(strcmp(time.hdr.label,'UADC016'))
        cfg = [];
        cfg.dataset = filename;
        cfg.continuous = 'yes';
        cfg.channel = 'UADC016'; % Time channel!
        light = ft_preprocessing(cfg);
        pix_sample = zeros(size(trig_sample.sample,1),1);
        ii = 0;
        for tt=1:length(light.trial)
            indt = find( diff(light.trial{tt}) < -1 );
            if ~isempty(indt)
                indt = indt([true, diff(indt) > light.fsample*0.15]);
                sampt = light.sampleinfo(tt,1) : light.sampleinfo(tt,2);
                pix_sample(ii+(1:nnz(indt)),1) = sampt(indt);
                ii = ii + nnz(indt);
            end
        end
        pix_sample = pix_sample(1:ii);
        if ~isempty(pix_sample)
            for tt = 1:length(pix_sample)
                d = pix_sample(tt) - trig_sample.sample ;
                d(d< 0) = 1000;
                [md, iid]  = min(d);
                if md < (light.fsample*0.05) % expect < 50ms delay
                    trig_sample.sample(iid) = pix_sample(tt);
                end
            end
        else
            trig_sample.sample = trig_sample.sample + round(0.018*light.fsample); % 18ms average delay
        end
    
    else
    
        trig_sample.sample = trig_sample.sample + round(0.018*time.fsample); % 18ms average delay
    end
end

%% Button presses
buttonpress = [];

if any(strcmp(time.hdr.label,'UADC006'))
    cfg = [];
    cfg.dataset = filename;
    cfg.continuous = 'yes';
    cfg.channel = {'UADC005';'UADC006';'UADC007';'UADC008'};
    buttons = ft_preprocessing(cfg);

    maxpeak = max(buttons.trial{1},[],2);

    for bb = 1:4

        if maxpeak(bb) > 2
            [PKS,LOCS] = findpeaks(buttons.trial{1}(bb,:),...
                'MinPeakDistance',0.01*buttons.fsample,...
                'MinPeakHeight',maxpeak(bb)*.7, 'MinPeakProminence',maxpeak(bb)*.7) ;
%             figure; plot(buttons.trial{1}(bb,:));
%             hold on; plot(LOCS,PKS,'r*')
%     
            buttonpress.(buttons.label{bb}) = LOCS;
        end
    end

    % old definition (edited on Nov 1 2024)
%     for bb = 1:4
%         ii = 0;
%         buttonpress.(buttons.label{bb}) = zeros(size(trig_sample,1),1);
% 
%         for tt=1:length(buttons.trial)
%             indt = find( diff(buttons.trial{tt}(bb,:),1,2) > 1 );
%             if ~isempty(indt)
%                 indt = indt([true, diff(indt) > buttons.fsample*0.2]);
%                 sampt = buttons.sampleinfo(tt,1) : buttons.sampleinfo(tt,2);
%                 buttonpress.(buttons.label{bb})(ii+(1:nnz(indt))) = sampt(indt);
%                 ii = ii + nnz(indt);
%             end
%         end
%         buttonpress.(buttons.label{bb}) = buttonpress.(buttons.label{bb})(1:ii);
%     end
end
end

