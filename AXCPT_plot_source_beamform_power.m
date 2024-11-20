clearvars
ftpath = '/home/liuzzil2/fieldtrip-20190812/';
addpath(ftpath)

ft_defaults
addpath ~/matlab_utils/

rootdir = '/data/EDB/AXCPT-Flanker3_Adult/';

subincl = readtable([rootdir,'/derivatives/meg_derivatives/MEG_data_inclusion.csv']);
sublist = subincl.SDAN( strcmp(subincl.AXCPT_usable,'yes' ) &  strcmp(subincl.coregistered,'yes' ) )  ;
n = length(sublist);

derdir = [rootdir, 'derivatives/meg_derivatives/'];

mu  =4;
gridres = 5; % 5mm grid

mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates
mri_mni.coordsys = 'mni';

load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

addpath('/data/liuzzil2/NIfTI/')
[m,ii]= min(sum(sourcemodel.pos.^2,2));
[x,y,z]=ind2sub(sourcemodel.dim,ii);

pow_A = zeros(nnz(sourcemodel.inside),n);
pow_B = zeros(nnz(sourcemodel.inside),n);
% pow_diff = zeros(nnz(sourcemodel.inside),n);
% pow_AX = zeros(nnz(sourcemodel.inside),n);
% pow_AY = zeros(nnz(sourcemodel.inside),n);
% pow_BX = zeros(nnz(sourcemodel.inside),n);
% pow_BY = zeros(nnz(sourcemodel.inside),n);

% pow_cue = zeros(nnz(sourcemodel.inside),n);
% pow_probe = zeros(nnz(sourcemodel.inside),n);


lowfv = [1, 4, 8, 13];
highv = [ 4, 8, 13, 30];

for ff = 1:4
   
lowf = lowfv(ff);
highf = highv(ff);

%%

for ss = 1:n
    fprintf('Loading subject %d/%d\n', ss, n)

    sub = num2str(sublist(ss));
    
    if strcmp(sub,'24531')
        processingfolder = [derdir,'/sub-',sub,'/ses-02/'];
    else
        processingfolder = [derdir,'/sub-',sub,'/'];
    end
    
    d = dir(processingfolder);
    d = struct2cell(d);
    filenames = d(1,3:end);
    
    df = find(contains(filenames,'task-axcpt_run-'));
    nrec = nnz(df);
    
    s_A = zeros(nnz(sourcemodel.inside),nrec);
    s_B = zeros(nnz(sourcemodel.inside),nrec);
%     s_diff = zeros(nnz(sourcemodel.inside),nrec);
%     s_AX = zeros(nnz(sourcemodel.inside),nrec);
%     s_AY = zeros(nnz(sourcemodel.inside),nrec);
%     s_BX = zeros(nnz(sourcemodel.inside),nrec);
%     s_BY = zeros(nnz(sourcemodel.inside),nrec);

%     s_cue = zeros(nnz(sourcemodel.inside),nrec);
%     s_probe = zeros(nnz(sourcemodel.inside),nrec);
    
    Ntrials = [];
    for ii = 1:nrec
       
        filename = filenames{df(ii)};
        datafile = sprintf('%s%s/Prep_Pow%d-%dHz_multiSpheres_%dmm_regmu%d.mat',...
                processingfolder,filename,lowf,highf,gridres,mu);
        if exist( datafile,'file')
            
            load(datafile)
       
            s_A(:,ii) = PcueA;
            s_B(:,ii) = PcueB;
%             s_diff(:,ii) = PdiffAB;
%             s_AX(:,ii) = PprobeAX;
%             s_AY(:,ii) = PprobeAY;
%             s_BX(:,ii) = PprobeBX;
%             s_BY(:,ii) = PprobeBY;
%             s_cue(:,ii) = Pcue;
%             s_probe(:,ii) = Pprobe;


            if ii == 1
                Ntrials = ntrials;
            else
                Ntrials.A   = cat(2, Ntrials.A, ntrials.A);
                Ntrials.B   = cat(2, Ntrials.B, ntrials.B);
%                 Ntrials.AX  = cat(2, Ntrials.AX, ntrials.AX);
%                 Ntrials.AY  = cat(2, Ntrials.AY, ntrials.AY);
%                 Ntrials.BX  = cat(2, Ntrials.BX, ntrials.BX);
%                 Ntrials.BY  = cat(2, Ntrials.BY, ntrials.BY);

%                 Ntrials.cues   = cat(2, Ntrials.cues, ntrials.cues);
%                 Ntrials.probes   = cat(2, Ntrials.probes, ntrials.probes);
            end
           
%            
        end
    end
    if ~isempty(Ntrials)
        % average over runs
        pow_A(:,ss)  = sum(s_A.*Ntrials.A ,2)/sum(Ntrials.A);
        pow_B(:,ss)  = sum(s_B.*Ntrials.B ,2)/sum(Ntrials.B);
%         pow_diff(:,ss)  = sum(s_diff.*Ntrials.B ,2)/sum(Ntrials.B);
%         pow_AX(:,ss) = sum(s_AX.*Ntrials.AX ,2)/sum(Ntrials.AX);
%         pow_AY(:,ss) = sum(s_AY.*Ntrials.AY ,2)/sum(Ntrials.AY);
%         pow_BX(:,ss) = sum(s_BX.*Ntrials.BX ,2)/sum(Ntrials.BX);
%         pow_BY(:,ss) = sum(s_BY.*Ntrials.BY ,2)/sum(Ntrials.BY);
    
%         pow_cue(:,ss)  = sum(s_cue.*Ntrials.cues ,2)/sum(Ntrials.cues);
%         pow_probe(:,ss) = sum(s_probe.*Ntrials.probes ,2)/sum(Ntrials.probes);
    end
    
end   

% nincl = sum(pow_A,1) ~= 0;
% pow_A = pow_A(:,nincl);
% pow_B = pow_B(:,nincl);
% pow_AX = pow_AX(:,nincl);
% pow_AY = pow_AY(:,nincl);
% pow_BX = pow_BX(:,nincl);
% pow_BY = pow_BY(:,nincl);
% 
% nincl = sum(pow_cue,1) ~= 0;
% pow_cue = pow_cue(:,nincl);
% pow_probe = pow_probe(:,nincl);

nincl = sum(pow_A,1) ~= 0;
pow_A = pow_A(:,nincl);
pow_B = pow_B(:,nincl);
% pow_diff = pow_diff(:,nincl);

%% Plot Freq

% for cc = 1:2
% 
% if cc ==1
    condname = 'AvsBcue';
%     x1 = pow_diff;

    x1 = pow_A;
    x2 = pow_B;
    name1 = 'Acue';
    name2 = 'Bcue';

% elseif cc ==2
%     condname = 'AXvsBXprobe';
%     x1 = pow_AX;
%     x2 = pow_BX;
%     name1 = 'AXprobe';
%     name2 = 'BXprobe';
% elseif cc ==3
%     condname = 'AXvsAYprobe';
%     x1 = pow_AX;
%     x2 = pow_AY;
%     name1 = 'AXprobe';
%     name2 = 'AYprobe';
% else
%     error('Unrecognized condition')
% end

% if cc ==1
%     condname = 'cues';
%     x1 = pow_cue;
%     
% elseif cc ==2
%     condname = 'probes';
%     x1 = pow_probe;
% end



% sp  = sqrt((std(x1,[],2).^2 + std(x2,[],2).^2) /2 );
% t = (mean(x1,2) - mean(x2,2)) ./ (sp*sqrt(2/n));

[~,~,~,stats]=ttest(x1',x2');
t = stats.tstat;
% 
% 
% T = zeros(sourcemodel.dim);
% T(sourcemodel.inside) = t;
% P = zeros(sourcemodel.dim);
% P(sourcemodel.inside) = mean(x1 - x2,2);
% 
% 
% 
% PT =cat(4,P,T);
% nii4d =  make_nii(PT, [1 1 1]*gridres ,[x,y,z]);
% save_nii(nii4d, sprintf('%s/grouplevel/AXCPT_source/%s_%d-%dHz.nii',derdir,condname,lowf,highf))

% [~,~,~,stats]=ttest(x1');
% t = stats.tstat;
T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = t;

P = zeros(sourcemodel.dim);
P(sourcemodel.inside) = mean(x1,2) - mean(x2,2) ;

sourceant =[];
sourceant.pow = P;
sourceant.T = T;

sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;

cfg = [];
cfg.parameter = {'pow','T'};
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.T(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';
    
%     sourceant.dim = gridl.dim;
%     sourceant.inside = gridl.inside;
%     sourceant.pos = gridl.pos;
%     cfg = [];
%     cfg.parameter = 'pow';
%     sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
%     sourceout_Int.pow(~sourceout_Int.inside) = 0;

% crang = [thresh max(sourceant.pow)];
crang = [];

mriplot = 'mni';
plotopt = 'ortho';  % 'slice'; %'ortho'


cfg = [];
cfg.method        = plotopt; %'ortho'
cfg.location   = 'min';

  
cfg.funparameter = 'pow';

cfg.maskparameter = 'T';
   
cfg.maskstyle     ='colormix';
%  cfg.opacitymap    ='vdown';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = [];
if strcmp(mriplot,'mni')
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
end

ft_sourceplot(cfg, sourceout_Int);

set(gcf,'position',[704   117   894   675],'color','w')
title(sprintf('%s %d-%dHz',condname,lowf,highf))


end

