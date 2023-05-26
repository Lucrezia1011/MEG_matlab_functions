
addpath /home/liuzzil2/fieldtrip-20190812/

ft_defaults
addpath ~/matlab_utils/

gridres = 5;
mu = 4;
mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates
    mri_mni.coordsys = 'mni';
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    sourcemodel.coordsys = 'mni';

subs = {'24531';'24563';'24590';'24580';'24626';'24581';'24482';...
    '24640';'24592';'24667';'24678' };


Theta_A = zeros(nnz(sourcemodel.inside),length(subs));
Theta_B = zeros(nnz(sourcemodel.inside),length(subs));
Theta_AX = zeros(nnz(sourcemodel.inside),length(subs));
Theta_AY = zeros(nnz(sourcemodel.inside),length(subs));
Theta_BX = zeros(nnz(sourcemodel.inside),length(subs));
Theta_BY = zeros(nnz(sourcemodel.inside),length(subs));

lowf = 4;
highf = 8;
%%

for ss = 1:length(subs)
    sub = subs{ss};
    if strcmp(sub,'24531')
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/ses-02/'];
    else
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
    end
    
    d = dir(processingfolder);
    d = struct2cell(d);
    filenames = d(1,3:end);
    
    df = find(contains(filenames,'task-axcpt_run-'));
    nrec = nnz(df);
    
    sERF_A = zeros(nnz(sourcemodel.inside),nrec);
    sERF_B = zeros(nnz(sourcemodel.inside),nrec);
    sERF_AX = zeros(nnz(sourcemodel.inside),nrec);
    sERF_AY = zeros(nnz(sourcemodel.inside),nrec);
    sERF_BX = zeros(nnz(sourcemodel.inside),nrec);
    sERF_BY = zeros(nnz(sourcemodel.inside),nrec);
    
    Ntrials = [];
    for ii = 1:nrec
       
        filename = filenames{df(ii)};
        datafile = sprintf('%s%s/Pow%d-%dHz_multiSpheres_%dmm_regmu%d.mat',...
            processingfolder,filename,lowf,highf,gridres,mu);
        if exist( datafile,'file')
            
            load(datafile)
       
            sERF_A(:,ii) = PcueA;
            sERF_B(:,ii) = PcueB;
            sERF_AX(:,ii) = PprobeAX;
            sERF_AY(:,ii) = PprobeAY;
            sERF_BX(:,ii) = PprobeBX;
            sERF_BY(:,ii) = PprobeBY;
%            
        end
    end
    

    Theta_A(:,ss)  = mean(sERF_A,2);
    Theta_B(:,ss)  = mean(sERF_B,2);
    Theta_AX(:,ss) = mean(sERF_AX,2);
    Theta_AY(:,ss) = mean(sERF_AY,2);
    Theta_BX(:,ss) = mean(sERF_BX,2);
    Theta_BY(:,ss) = mean(sERF_BY,2);
    
    
end   


%% Plot Freq

crang = [];

mriplot = 'mni';
plotopt = 'slice';  % 'slice'; %'ortho'

T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = mean(Theta_AY,2);% -  mean(Theta_AX,2) ; 
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