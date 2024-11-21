function combineDs(filenames, basenfile, datasetname, datapath)

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
cd(datapath)

ds =  cell(1,length(filenames));
data = [];
ntrials = 0;
fileSize = 0;
for ii = 1:length(filenames)
    ds{ii}=readCTFds(filenames{ii});
    data_raw = ft_read_data(filenames{ii});
    data = cat(3, data, data_raw);
    ntrials = ntrials + ds{ii}.res4.no_trials;
    fileSize = fileSize +  ds{ii}.meg4.fileSize;
end
data = permute(data, [2,1,3]);
ds_out = ds{basenfile};

ds_out.res4.no_trials = ntrials;
ds_out.meg4.fileSize = fileSize;


for m = 1:length(ds_out.mrk)
    mrk_name = ds_out.mrk(m).Name;
    
    ds_out.mrk(m).trial = [];
    ds_out.mrk(m).time = [];
    mtrials = 0;
    for ii = 1:length(filenames)
        indf = 0;
        for ms = 1:length(ds{ii}.mrk)
            if strcmp(ds{ii}.mrk(ms).Name, mrk_name)
                indf = ms;
            end
        end
        ds_out.mrk(m).trial = cat(2, ds_out.mrk(m).trial, ds{ii}.mrk(indf).trial + mtrials);
        ds_out.mrk(m).time = cat(2, ds_out.mrk(m).time, ds{ii}.mrk(indf).time);
        
        mtrials = mtrials + ds{ii}.res4.no_trials;
    end
    
end

dse=writeCTFds(datasetname,ds_out,data,'T');
writeMarkerFile(MarkerFile,marker)
