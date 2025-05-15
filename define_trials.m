function [datanew,ttdel] = define_trials(samples,data,trialed,opts)
% Lucrezia Liuzzi, last updated 2024/04/05
% 
% Defines trials compatible with field_trip
% [datanew,ttdel] = define_trials(samples,data,trialed,opts)
% samples   = vector with data samples to be new t=0  
% data      = continuos data in fieldtrip format
% samples, % NO TIME INPUT NEEDED, EDIT ON APRIL 5h 2024
%             used to identify discontinuities (output of matchTriggers)
% trialed   = new trial time edges in seconds, e.g. [-2 2]
% opts      = verbose (true,false)
l = length(samples);

f = data.fsample;
wind_press = trialed(1)*f+1:trialed(2)*f; 

n = length(wind_press);
datanew = data;

datanew.trial = cell(1,l);
datanew.time = cell(1,l);
datanew.sampleinfo = repmat([1 n],l,1);
% datanew.sampleinfo = zeros(l,2);
time  = data.time{1};
tjump = find(diff(time)>=(2/data.fsample)); 
% tjump = find(abs(diff(time))>1); % Old method, using samples from matchTriggers, edit on April 5th 2024
tjump(end+1) = length(time);

ttdel = [];
if opts == true
    fprintf('Defining %d trials between %.2fs and %.2fs of samples\n',l,trialed(1),trialed(2))
end
for tt= 1:l
    sampwind = round(samples(tt)+wind_press);
    if sampwind(1)<0 || sampwind(end)>data.sampleinfo(2)
        ttdel = cat(1,ttdel,tt);
        warning('Requesting time window outside data bounds\nSkipping trial %d\n',tt)
    elseif isempty(intersect(tjump,sampwind))
        datanew.trial{tt} = data.trial{1}(:,sampwind);
        datanew.time{tt}  = linspace(wind_press(1)/f,wind_press(end)/f, n);
%         datanew.sampleinfo(tt,:)  = [sampwind(1) sampwind(end)];
    else
        ttdel = cat(1,ttdel,tt);
        warning('Time discontinuity found in trial %d\n',tt)
    end
    
end

datanew.trial(ttdel) =[];
datanew.time(ttdel) = [];
datanew.sampleinfo(ttdel,:) = [];
if opts ==true
    fprintf('New data contains %d trials\n',l-length(ttdel))
end



