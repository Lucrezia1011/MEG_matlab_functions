
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

subs = {'24531';'24563';'24590';'24580';'24588';'24626';'24581';'24482';...
    '24640';'24592';'24667';'24678';'24663' };


Theta_Flan = zeros(nnz(sourcemodel.inside),length(subs)*2);
Theta_Flan_C = zeros(nnz(sourcemodel.inside),length(subs)*2);
Theta_Flan_I = zeros(nnz(sourcemodel.inside),length(subs)*2);
Theta_Flan_err = zeros(nnz(sourcemodel.inside),length(subs)*2);

Alpha_Flan = zeros(nnz(sourcemodel.inside),length(subs)*2);
Alpha_Flan_C = zeros(nnz(sourcemodel.inside),length(subs)*2);
Alpha_Flan_I = zeros(nnz(sourcemodel.inside),length(subs)*2);

ERF_err = zeros(nnz(sourcemodel.inside),780,length(subs)*2);
ERF_Icorr = zeros(nnz(sourcemodel.inside),780,length(subs)*2);
ERF_Ccorr = zeros(nnz(sourcemodel.inside),780,length(subs)*2);
ERF_I = zeros(nnz(sourcemodel.inside),600,length(subs)*2);
ERF_C = zeros(nnz(sourcemodel.inside),600,length(subs)*2);

n = 0;
for ss = 1:length(subs)
    sub = subs{ss};
    datapath = ['/data/EDB/MEG_AXCPT_Flanker/data/sub-',sub,'/meg/'];

    if strcmp(sub,'24531')
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/ses-02/'];
    else
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
    end
    for ii = 1:2
        if strcmp(sub,'24531')
            filename = ['sub-',sub,'_ses-02_task-flanker_run-',num2str(ii),'_meg.ds'];
        else
            filename = ['sub-',sub,'_task-flanker_run-',num2str(ii),'_meg.ds'];
        end
       
        n = n+1;

        load( sprintf('%s%s/Congruency_multiSpheres_%dmm_regmu%d.mat',...
            processingfolder,filename(1:end-3),gridres,mu) )
    
        Theta_Flan(:,n) = Theta_Flan_Inc_Cong;
        Theta_Flan_C(:,n) = Theta_Flan_Cong;
        Theta_Flan_I(:,n) = Theta_Flan_Inc;
        Theta_Flan_err(:,n) = Theta_Resp_Err_Corr;
        
        Alpha_Flan(:,n) = Alpha_Flan_Inc_Cong;
        Alpha_Flan_C(:,n) = Alpha_Flan_Cong;
        Alpha_Flan_I(:,n) = Alpha_Flan_Inc;

        load( sprintf('%s%s/ERN_multiSpheres_%dmm_regmu%d.mat',...
            processingfolder,filename(1:end-3),gridres,mu));
        
        ERF_err(:,:,n) = (Pcomm);
        ERF_Icorr(:,:,n) = (PIcorr);
        ERF_Ccorr(:,:,n) = (PCcorr);
        
        ERF_I(:,:,n) = (PIflan);
        ERF_C(:,:,n) = (PCflan);
        
    end
end   


%% Sign flip
nv = size(ERF_err,1);

iitb = resptime >= -0.5 & resptime <= -0.2;
ERF_err = ERF_err - mean(ERF_err(:,iitb,:),2);
ERF_Ccorr = ERF_Ccorr - mean(ERF_Ccorr(:,iitb,:),2);
ERF_Icorr = ERF_Icorr - mean(ERF_Icorr(:,iitb,:),2);

iitb = flantime  <= 0;
ERF_I = ERF_I - mean(ERF_I(:,iitb,:),2);
ERF_C = ERF_C - mean(ERF_C(:,iitb,:),2);
for ix = 1:nv
    erfcon = [squeeze(ERF_I(ix,:,:)); squeeze(ERF_C(ix,:,:))];
    vec = corr(erfcon);
    ERF_I(ix,:,:) = squeeze(ERF_I(ix,:,:)).*sign(vec(1,:));
    ERF_C(ix,:,:) = squeeze(ERF_C(ix,:,:)).*sign(vec(1,:));
    
    
    erfcon = [squeeze(ERF_Ccorr(ix,:,:)); squeeze(ERF_Icorr(ix,:,:)) ; squeeze(ERF_err(ix,:,:))];
    vec = corr(erfcon);

    ERF_err(ix,:,:) = squeeze(ERF_err(ix,:,:)).*sign(vec(1,:));
    ERF_Icorr(ix,:,:) = squeeze(ERF_Icorr(ix,:,:)).*sign(vec(1,:));
    ERF_Ccorr(ix,:,:) = squeeze(ERF_Ccorr(ix,:,:)).*sign(vec(1,:));
   
end


%% ROI analysis
atlas =  ft_read_mri('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates
mri_mni.coordsys = 'mni';

%% Plot ERN

close all
% t0 = 0.05; t1 = 0.12;
t0 = 0.25; t1 = 0.5;

    iit = resptime >= t0 & resptime <= t1;
%        iit = flantime >= t0 & flantime <= t1;
       
    mriplot = 'mni';
    plotopt = 'ortho';  % 'slice'; %'ortho'

%     ERF = (ERF_Icorr + ERF_Ccorr)/2;
%     xcorr = abs( mean(ERF(:,flantime >= -0.5 & flantime <= 0,:),2) );
%     xerr = abs( mean( ERF(:,iit,:),2));

    iitb = resptime >= -0.5 & resptime <= -0.25;
%     iitb = flantime <= 0;
    x1 = ( mean( ERF_err(:,iit,:),2));
    x2 = ( mean(ERF_Icorr(:,iit,:),2) );
%     x2 = ( mean(ERF(:,iitb,:),2) );

    sp  = sqrt((std(x1,[],3).^2 + std(x2,[],3).^2) /2 );
%     t = (abs(mean(x1,3)) - abs(mean(x2,3))) ./ (sp*sqrt(2/n));
    t = abs(mean(x1,3) - mean(x2,3)) ./ (sp*sqrt(2/n));
    
    
    

    T = zeros(sourcemodel.dim);
    T(sourcemodel.inside) = t;% (abs(mean(xerr,3)))- abs(mean(xcorr,3))) * 1e15;
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

    sourcet = sourceant;
    T = zeros(sourcemodel.dim);
    T(sourcemodel.inside) = abs(t);
    sourcet.pow = T;

    sourcet_Int  = ft_sourceinterpolate(cfg, sourcet , mri_mni);

    sourceout_Int.T = abs(sourcet_Int.pow ) ;
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

    crang = [];
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'T';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = [];
    if strcmp(mriplot,'mni')
        cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
    end

    ft_sourceplot(cfg, sourceout_Int);
    set(gcf,'position',[704   117   894   675],'color','w')
    title(sprintf('error - correct %d-%d ms',round(t0*1000),round(t1*1000)))
%     title(sprintf('flanker - baseline %d-%d ms',round(t0*1000),round(t1*1000)))

%%
mnicor =  [-54 -31 53];
[xdiff,iix] = min(sqrt( sum((sourceant.pos - mnicor/10).^2 ,2) ) );

T = zeros(sourcemodel.dim);
T(iix) = 1;
T = T(sourcemodel.inside);
ix = find(T);

erfcon = [squeeze(ERF_I(ix,:,:)); squeeze(ERF_C(ix,:,:))];
vec = corr(erfcon);
% erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule
s = 20;
figure; set(gcf,'color','w','position', [83   139   1575  433])

subplot(121); co = get(gca, 'colororder');
plot(flantime, smooth(mean(squeeze(ERF_C(ix,:,:)).*sign(vec(1,:)),2),s),'linewidth',2 )
hold all
plot(flantime, smooth(mean(squeeze(ERF_I(ix,:,:)).*sign(vec(1,:)),2),s),'linewidth',2 )
fill([flantime, fliplr(flantime)] , ...
    [smooth(mean(squeeze(ERF_C(ix,:,:)).*sign(vec(1,:)),2),s) + smooth(std(squeeze(ERF_C(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s); ...
    flipud(smooth(mean(squeeze(ERF_C(ix,:,:)).*sign(vec(1,:)),2),s) - smooth(std(squeeze(ERF_C(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s))], ...
    co(1,:),'edgecolor','none','facealpha',0.3)
fill([flantime, fliplr(flantime)] , ...
    [smooth(mean(squeeze(ERF_I(ix,:,:)).*sign(vec(1,:)),2),s) + smooth(std(squeeze(ERF_I(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s); ...
    flipud(smooth(mean(squeeze(ERF_I(ix,:,:)).*sign(vec(1,:)),2),s) - smooth(std(squeeze(ERF_I(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s))], ...
    co(2,:),'edgecolor','none','facealpha',0.3)
grid on
% ylim([-1 1]*1e-13)


erfcon = [squeeze(ERF_Icorr(ix,:,:)); squeeze(ERF_Ccorr(ix,:,:)) ; squeeze(ERF_err(ix,:,:))];
vec = corr(erfcon);
% erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule
legend('cong','incng','location','best')
title(sprintf('Flanker, mni[%d, %d, %d]',mnicor(1),mnicor(2),mnicor(3)))

 subplot(122)
plot(resptime, smooth(mean(squeeze(ERF_Ccorr(ix,:,:)).*sign(vec(1,:)),2),s),'linewidth',2 )
hold all
plot(resptime,smooth( mean(squeeze(ERF_Icorr(ix,:,:)).*sign(vec(1,:)),2),s),'linewidth',2 )
plot(resptime, smooth(mean(squeeze(ERF_err(ix,:,:)).*sign(vec(1,:)),2),s) ,'k' )

fill([resptime, fliplr(resptime)] , ...
    [smooth(mean(squeeze(ERF_Ccorr(ix,:,:)).*sign(vec(1,:)),2),s) + smooth(std(squeeze(ERF_Ccorr(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s); ...
    flipud(smooth(mean(squeeze(ERF_Ccorr(ix,:,:)).*sign(vec(1,:)),2),s) - smooth(std(squeeze(ERF_Ccorr(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s))], ...
    co(1,:),'edgecolor','none','facealpha',0.3)
fill([resptime, fliplr(resptime)] , ...
    [smooth(mean(squeeze(ERF_Icorr(ix,:,:)).*sign(vec(1,:)),2),s) + smooth(std(squeeze(ERF_Icorr(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s); ...
    flipud(smooth(mean(squeeze(ERF_Icorr(ix,:,:)).*sign(vec(1,:)),2),s) - smooth(std(squeeze(ERF_Icorr(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s))], ...
    co(2,:),'edgecolor','none','facealpha',0.3)
fill([resptime, fliplr(resptime)] , ...
    [smooth(mean(squeeze(ERF_err(ix,:,:)).*sign(vec(1,:)),2),s) + smooth(std(squeeze(ERF_err(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s); ...
    flipud(smooth(mean(squeeze(ERF_err(ix,:,:)).*sign(vec(1,:)),2),s) - smooth(std(squeeze(ERF_err(ix,:,:)).*sign(vec(1,:)),[],2)/sqrt(n),s))], ...
    [0.5 0.5 0.5],'edgecolor','none','facealpha',0.3)


grid on
% ylim([-1 1]*1e-13)

legend('cong','incng','err','location','best')
% set(gcf,'name',sprintf('Corr ERF peak at %dms ',round(tt*1e3)))
% saveas(gcf,sprintf('%s/Comm_ERF_%dms_%s.jpg',resultsfolder,round(tt*1e3),stimname))

% title(sprintf('Incong Comm %d,%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm_%d-%dms_resp.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000))
% 
% saveas(gcf,sprintf('%s/Incong_Corr-base_%d-%dms_13-30Hz_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000,stimname))

title(sprintf('Response, mni[%d, %d, %d]',mnicor(1),mnicor(2),mnicor(3)))

 %% Plot motor response
iit = resptime >= -0.05 & resptime <= 0.1;
ERFresp = (ERF_Icorr + ERF_Ccorr)/2;

xmot = abs( mean(  ERFresp(:,iit,:)  ,2) );

iit = resptime >= -0.5 & resptime <= -0.35;
xbase = abs( mean( ERFresp(:,iit,:),2));

sp  = sqrt((std(xmot,[],3).^2 + std(xbase,[],3).^2) /2 );
t = (mean(xmot,3) - mean(xbase,3)) ./ (sp*sqrt(2/n)); 

crang = [];
%%
mriplot = 'mni';
plotopt = 'ortho';  % 'slice'; %'ortho'

T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = (mean(xmot,3) - mean(xbase,3)) * 1e15; 
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
    
sourcet = sourceant;
T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = abs(t); 
sourcet.pow = T;

sourcet_Int  = ft_sourceinterpolate(cfg, sourcet , mri_mni);

sourceout_Int.T = abs(sourcet_Int.pow ) ;
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
cfg.maskparameter = 'T';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = [3 8];
if strcmp(mriplot,'mni')
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
end

ft_sourceplot(cfg, sourceout_Int);
colormap default

set(gcf,'position',[704   117   894   675],'color','w')
%%
mnicor =  [-37 -25 55]/10;
[xdiff,iix] = min(sqrt( sum((sourceant.pos - mnicor).^2 ,2) ) );

T = zeros(sourcemodel.dim);
T(iix) = 1;
T = T(sourcemodel.inside);
ix = find(T);

erfcon = squeeze(ERFresp(ix,:,:));
vec = corr(erfcon);
erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule

figure; hold all
for ii = 1:length(subs)
    plot(resptime, mean(erfcon(:,(ii-1)*2 + (1:2)),2) )
end

 %% Plot Flanker stimulus
flantime = resptime(181:end);
 
iit = flantime >= 0.2 & flantime <= 0.4;
ERFflan = (ERF_I + ERF_C)/2;

xvis = abs( mean(  ERFflan(:,iit,:)  ,2) );

iit = flantime >= -0.2 & flantime <= -0.0;
xbase = abs( mean( ERFflan(:,iit,:),2));

sp  = sqrt((std(xvis,[],3).^2 + std(xbase,[],3).^2) /2 );
t = (mean(xvis,3) - mean(xbase,3)) ./ (sp*sqrt(2/n)); 

crang = [];
%%
mriplot = 'mni';
plotopt = 'ortho';  % 'slice'; %'ortho'

T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = (mean(xvis,3) - mean(xbase,3)) * 1e15; 
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
    
sourcet = sourceant;
T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = abs(t); 
sourcet.pow = T;

sourcet_Int  = ft_sourceinterpolate(cfg, sourcet , mri_mni);

sourceout_Int.T = abs(sourcet_Int.pow ) ;
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
cfg.maskparameter = 'T';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = [0 4];
if strcmp(mriplot,'mni')
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
end

ft_sourceplot(cfg, sourceout_Int);
colormap default

set(gcf,'position',[704   117   894   675],'color','w')
%%
mnicor =  [2 41 24]/10;
[xdiff,iix] = min(sqrt( sum((sourceant.pos - mnicor).^2 ,2) ) );

T = zeros(sourcemodel.dim);
T(iix) = 1;
T = T(sourcemodel.inside);
ix = find(T);

erfcon = squeeze(ERFflan(ix,:,:));
vec = corr(erfcon);
erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule

figure;
plot(flantime, squeeze(ERFflan(ix,:,:)).*sign(vec(1,:)) )


%% Plot Freq

% Theta_Flan(:,n) = Theta_Flan_Inc_Cong;
%         Theta_Flan_C(:,n) = Theta_Flan_Cong;
%         Theta_Flan_I(:,n) = Theta_Flan_Inc;
%         Theta_Flan_err(:,n) = Theta_Resp_Err_Corr;
%         
%         Alpha_Flan(:,n) = Alpha_Flan_Inc_Cong;
%         Alpha_Flan_C(:,n) = Alpha_Flan_Cong;
%         Alpha_Flan_I(:,n) = Alpha_Flan_Inc;

crang = [];

mriplot = 'mni';
plotopt = 'ortho';  % 'slice'; %'ortho'

T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = mean(Alpha_Flan ,2); 
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
cfg.location   = 'min';

cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
if strcmp(mriplot,'mni')
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
end

ft_sourceplot(cfg, sourceout_Int);
colormap default

set(gcf,'position',[704   117   894   675],'color','w')