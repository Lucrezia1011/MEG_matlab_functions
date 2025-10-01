ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
addpath(ftpath)
ft_defaults

roothpath = '/data/EDB/MEG_AXCPT_Flanker/';

addpath ~/matlab_utils/

gridres = 5;
mu = 4;
mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates
mri_mni.coordsys = 'mni';

load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

sublist = {'24531';'24563';'24590';'24580';'24626';'24581';'24482';...
    '24640';'24592';'24667';'24678' };


pow_A = zeros(nnz(sourcemodel.inside),length(sublist));
pow_B = zeros(nnz(sourcemodel.inside),length(sublist));
pow_AX = zeros(nnz(sourcemodel.inside),length(sublist));
pow_AY = zeros(nnz(sourcemodel.inside),length(sublist));
pow_BX = zeros(nnz(sourcemodel.inside),length(sublist));
pow_BY = zeros(nnz(sourcemodel.inside),length(sublist));

lowf = 4;
highf = 8;
%%

for ss = 1:length(sublist)
    sub = sublist{ss};
    if strcmp(sub,'24531')
        processingfolder = [roothpath,'derivatives/sub-',sub,'/ses-02/'];
    else
        processingfolder = [roothpath,'derivatives/sub-',sub,'/'];
    end
    
    d = dir(processingfolder);
    d = struct2cell(d);
    filenames = d(1,3:end);
    
    df = find(contains(filenames,'task-axcpt_run-'));
    nrec = nnz(df);
    
    s_A = zeros(nnz(sourcemodel.inside),nrec);
    s_B = zeros(nnz(sourcemodel.inside),nrec);
    s_AX = zeros(nnz(sourcemodel.inside),nrec);
    s_AY = zeros(nnz(sourcemodel.inside),nrec);
    s_BX = zeros(nnz(sourcemodel.inside),nrec);
    s_BY = zeros(nnz(sourcemodel.inside),nrec);
    
    Ntrials = [];
    for ii = 1:nrec
       
        filename = filenames{df(ii)};
        datafile = sprintf('%s%s/Pow%d-%dHz_multiSpheres_%dmm_regmu%d.mat',...
            processingfolder,filename,lowf,highf,gridres,mu);
        if exist( datafile,'file')
            
            load(datafile)
       
            s_A(:,ii) = PcueA;
            s_B(:,ii) = PcueB;
            s_AX(:,ii) = PprobeAX;
            s_AY(:,ii) = PprobeAY;
            s_BX(:,ii) = PprobeBX;
            s_BY(:,ii) = PprobeBY;
%            
        end
    end
    
    % average over runs
    pow_A(:,ss)  = mean(s_A,2);
    pow_B(:,ss)  = mean(s_B,2);
    pow_AX(:,ss) = mean(s_AX,2);
    pow_AY(:,ss) = mean(s_AY,2);
    pow_BX(:,ss) = mean(s_BX,2);
    pow_BY(:,ss) = mean(s_BY,2);
    
    
end   


%% Plot Freq

% condition to plot
P = mean(pow_AY,2) -  mean(pow_AX,2) ; 


crang = [];

mriplot = 'mni';
plotopt = 'slice';  % 'slice'; %'ortho'

T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = P;
sourceant =[];
sourceant.pow = T;
    
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;

cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';
    
%     sourceant.dim = gridl.dim;
%     sourceant.inside = gridl.inside;
%     sourceant.pos = gridl.pos;
%     cfg = [];
%     cfg.parameter = 'pow';
%     sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
%     sourceout_Int.pow(~sourceout_Int.inside) = 0;
    


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
