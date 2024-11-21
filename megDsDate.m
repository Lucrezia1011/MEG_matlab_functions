function TaskDate = megDsDate(datapath,filename)
% Get collection date of MEG dataset
cd([datapath,filename])

fileID = fopen([filename(1:end-2),'infods'],'r');
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