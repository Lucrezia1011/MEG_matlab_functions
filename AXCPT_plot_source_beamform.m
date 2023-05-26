
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

% Theta_A = zeros(nnz(sourcemodel.inside),length(subs)*2);
% Theta_B = zeros(nnz(sourcemodel.inside),length(subs)*2);
% Theta_AX = zeros(nnz(sourcemodel.inside),length(subs)*2);
% Theta_AY = zeros(nnz(sourcemodel.inside),length(subs)*2);
% Theta_BX = zeros(nnz(sourcemodel.inside),length(subs)*2);
% Theta_BY = zeros(nnz(sourcemodel.inside),length(subs)*2);

ERF_A  = zeros(nnz(sourcemodel.inside),900,length(subs));
ERF_B  = zeros(nnz(sourcemodel.inside),900,length(subs));
ERF_AX = zeros(nnz(sourcemodel.inside),900,length(subs));
ERF_AY = zeros(nnz(sourcemodel.inside),900,length(subs));
ERF_BX = zeros(nnz(sourcemodel.inside),900,length(subs));
ERF_BY = zeros(nnz(sourcemodel.inside),900,length(subs));
%%

for ss = 1:length(subs)
    sub = subs{ss};
    if ss == 1
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/ses-02/'];
    else
        processingfolder = ['/data/EDB/MEG_AXCPT_Flanker/derivatives/sub-',sub,'/'];
    end
    
    d = dir(processingfolder);
    d = struct2cell(d);
    filenames = d(1,3:end);
    
    df = find(contains(filenames,'task-axcpt_run-'));
    nrec = nnz(df);
    
    sERF_A = zeros(nnz(sourcemodel.inside),900,nrec);
    sERF_B = zeros(nnz(sourcemodel.inside),900,nrec);
    sERF_AX = zeros(nnz(sourcemodel.inside),900,nrec);
    sERF_AY = zeros(nnz(sourcemodel.inside),900,nrec);
    sERF_BX = zeros(nnz(sourcemodel.inside),900,nrec);
    sERF_BY = zeros(nnz(sourcemodel.inside),900,nrec);
    
    Ntrials = [];
    for ii = 1:nrec
       
        filename = filenames{df(ii)};
       
        if exist( sprintf('%s%s/ERF_multiSpheres_%dmm_regmu%d.mat',...
            processingfolder,filename,gridres,mu),'file')
            
            load( sprintf('%s%s/ERF_multiSpheres_%dmm_regmu%d.mat',...
                processingfolder,filename,gridres,mu))

            
            sERF_A(:,:,ii)= PerpA;
            sERF_B(:,:,ii) = PerpB;
            sERF_AX(:,:,ii) = PerpAX;
            sERF_AY(:,:,ii) = PerpAY;
            sERF_BX(:,:,ii) = PerpBX;
            sERF_BY(:,:,ii) = PerpBY;
            
            if ii == 1
                Ntrials = ntrials;
            else
                Ntrials.A   = cat(2, Ntrials.A, ntrials.A);
                Ntrials.B   = cat(2, Ntrials.B, ntrials.B);
                Ntrials.AX  = cat(2, Ntrials.AX, ntrials.AX);
                Ntrials.AY  = cat(2, Ntrials.AY, ntrials.AY);
                Ntrials.BX  = cat(2, Ntrials.BX, ntrials.BX);
                Ntrials.BY  = cat(2, Ntrials.BY, ntrials.BY);
            end
% 
%             load( sprintf('%s%s/Theta_multiSpheres_%dmm_regmu%d.mat',...
%                 processingfolder,filename(1:end-3),gridres,mu))
% 
%             Theta_A(:,n) =PcueA;
%             Theta_B(:,n) = PcueB;
%             Theta_AX(:,n) =PprobeAX;
%             Theta_AY(:,n) = PprobeAY;
%             Theta_BX(:,n) = PprobeBX;
%             Theta_BY(:,n) = PprobeBY;
%            
        end
    end
    
    
    svec = zeros(size(sERF_A,1),nrec);
    [~,im]= max(Ntrials.AX);
    for ix = 1:size(sERF_A,1)
        ca = corr(squeeze(sERF_A(ix,:,:)));
        cb = corr(squeeze(sERF_B(ix,:,:)));
        
        cx = corr(squeeze(sERF_AX(ix,:,:)));
        cy = corr(squeeze(sERF_AY(ix,:,:)));
%         corr(squeeze(sERF_BX(ix,:,:))) % fewer trials
%         corr(squeeze(sERF_BY(ix,:,:)))

        svec(ix,:) = sign(mean([ca(im,:);cb(im,:);cx(im,:);cy(im,:)],1));
    end
    
    svec = reshape(svec,[size(sERF_A,1),1,nrec]);
    
    
%     sERF_A1  = sERF_AX.*reshape(Ntrials.A,[1,1,nrec]);
%     figure; plot(erptime,squeeze(sERF_A1(100,:,:)))
%   
%     sERF_A2 = sERF_AX.*reshape(Ntrials.A,[1,1,nrec]).*svec;
%     figure; plot(erptime,squeeze(sERF_A2(100,:,:)))
%     
%     sERF_A3 = sum(sERF_AX.*reshape(Ntrials.A,[1,1,nrec]).*svec ,3);
%     hold ; plot(erptime,squeeze(sERF_A3(100,:,:)),'k','linewidth',2)
%     
%     sERF_A4 = sum(sERF_AX.*reshape(Ntrials.A,[1,1,nrec]).*svec ,3) / sum(Ntrials.A);
%     hold on; plot(erptime,squeeze(sERF_A4(100,:)),'k','linewidth',2)
%     
    mERF_A  = sum(sERF_A.*reshape(Ntrials.A,[1,1,nrec]).*svec ,3) / sum(Ntrials.A);
    mERF_B  = sum(sERF_B.*reshape(Ntrials.B,[1,1,nrec]).*svec ,3) / sum(Ntrials.B);
    mERF_AX = sum(sERF_AX.*reshape(Ntrials.AX,[1,1,nrec]).*svec ,3) / sum(Ntrials.AX);
    mERF_AY = sum(sERF_AY.*reshape(Ntrials.AY,[1,1,nrec]).*svec ,3) / sum(Ntrials.AY);
    mERF_BX = sum(sERF_BX.*reshape(Ntrials.BX,[1,1,nrec]).*svec ,3) / sum(Ntrials.BX);
    mERF_BY = sum(sERF_BY.*reshape(Ntrials.BY,[1,1,nrec]).*svec ,3) / sum(Ntrials.BY);

%     figure(ss);clf; 
%     subplot(121); plot(erptime,squeeze(sERF_A(100,:,:)))
%     hold on
%     subplot(122); plot(erptime,squeeze(sERF_AX(100,:,:)))
%     hold on
%     
%     subplot(121); plot(erptime,mERF_A(100,:),'k','linewidth',2)
%     title(sub)
%     subplot(122); plot(erptime,mERF_AX(100,:),'k','linewidth',2)
%     drawnow
    ERF_A(:,:,ss)  = mERF_A;
    ERF_B(:,:,ss)  = mERF_B;
    ERF_AX(:,:,ss) = mERF_AX;
    ERF_AY(:,:,ss) = mERF_AY;
    ERF_BX(:,:,ss) = mERF_BX;
    ERF_BY(:,:,ss) = mERF_BY;
    
    
end   

%% Sign flip
nv = size(ERF_A,1);
sERF_A = cell(nv,1);
sERF_B = cell(nv,1);
sERF_AX = cell(nv,1);
sERF_AY = cell(nv,1);
sERF_BX = cell(nv,1);
sERF_BY = cell(nv,1);

for ix = 1:nv
    erfcon = [squeeze(ERF_A(ix,:,:)); squeeze(ERF_B(ix,:,:))];
    vec = corr(erfcon);
    sERF_A{ix} = squeeze(ERF_A(ix,:,:)).*sign(vec(1,:));
    sERF_B{ix} = squeeze(ERF_B(ix,:,:)).*sign(vec(1,:));
    
    
    erfcon = [squeeze(ERF_AX(ix,:,:)); squeeze(ERF_AY(ix,:,:)) ; squeeze(ERF_BX(ix,:,:)); squeeze(ERF_BY(ix,:,:))];
    vec = corr(erfcon);

    sERF_AX{ix} = squeeze(ERF_AX(ix,:,:)).*sign(vec(1,:));
    sERF_AY{ix} = squeeze(ERF_AY(ix,:,:)).*sign(vec(1,:));
    sERF_BX{ix} = squeeze(ERF_BX(ix,:,:)).*sign(vec(1,:));
    sERF_BY{ix} = squeeze(ERF_BY(ix,:,:)).*sign(vec(1,:));
    
end

for ix = 1:nv
    ERF_AX(ix,:,:) = sERF_AX{ix};
    ERF_AY(ix,:,:) = sERF_AY{ix};
    ERF_BX(ix,:,:) = sERF_BX{ix};
    ERF_BY(ix,:,:) = sERF_BY{ix};
    
    ERF_A(ix,:,:) = sERF_A{ix};
    ERF_B(ix,:,:) = sERF_B{ix};


end
%% Plot ERF
close all
% for tx = 0:0.05:0.4
%     iit = erptime >= tx & erptime <= (tx+0.1);
    iit = erptime >= 0.25 & erptime <= 0.4;

    n = length(subs);
    
%     t = 0.4;
%     [~,iit] = min(abs(erptime - t));

     % contrast with two conditions
    x1 = ( mean(ERF_B(:,iit,:),2) );
    x2 = ( mean( ERF_A(:,iit,:),2));

    sp  = sqrt((std(x1,[],3).^2 + std(x2,[],3).^2) /2 );
    t = (abs(mean(x1,3)) - abs(mean(x2,3))) ./ (sp*sqrt(2/n));


    % contrast with baseline
%     iitb = erptime >= -0.5 & erptime <= -0.2;
%     x1 = ( mean(ERF_AY(:,iit,:),2) );
%     x2 = ( mean( ERF_AY(:,iitb,:),2));
% 
%     sp  = sqrt((std(x1,[],3).^2 + std(x2,[],3).^2) /2 );
%     t = (abs(mean(x1,3)) - abs(mean(x2,3))) ./ (sp*sqrt(2/n));
%     
    
%     iitb = erptime >= -0.5 & erptime <= -0.1;
    mriplot = 'mni';
 
% 
%     x1 = abs( mean(ERF_AY(:,iit,:),2) );
%     x2 = abs( mean( ERF_AX(:,iitb,:),2));
% %     x2 = abs( mean( ERF_B(:,iit,:),2));
%     
%     sp  = sqrt((std(x1,[],3).^2 + std(x2,[],3).^2) /2 );
%     t = (mean(x1,3) - mean(x2,3)) ./ (sp*sqrt(2/n));

    

    T = zeros(sourcemodel.dim);
    T(sourcemodel.inside) = t; %mean(abs(x1 - x2),3) * 1e15;
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



    crang = [-4 4];
    cfg = [];
    cfg.method        = 'ortho'; %'ortho'
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
    set(gcf,'position',[704   117   894   675],'color','w')

%     title(sprintf('B-A %d-%d ms',round(tx*1e3), round((tx+0.1)*1e3)))
    drawnow
% colormap default


%% ERF plot with sign
 iit = erptime >= -0.5 & erptime <= 1;

mnicor =  [-7 40 40]/10;
[xdiff,iix] = min(sqrt( sum((sourceant.pos - mnicor).^2 ,2) ) );
s =20;

T = zeros(sourcemodel.dim);
T(iix) = 1;
T = T(sourcemodel.inside);
ix = find(T);

erfcon = [squeeze(ERF_A(ix,iit,:)); squeeze(ERF_B(ix,iit,:))];
vec = corr(erfcon);
% erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule
figure; set(gcf,'color','w','position', [83   139   1575  433])
subplot(121)
A = squeeze(ERF_A(ix,:,:)).*sign(vec(1,:));
B = squeeze(ERF_B(ix,:,:)).*sign(vec(1,:));

plot(erptime, smooth(mean(A,2),s),'linewidth',2 )
hold all
plot(erptime,smooth( mean(B,2),s),'linewidth',2 )
grid on
c = get(gca,'colororder');
fill([erptime, fliplr(erptime)],...
    [smooth(mean(A,2) + std(A,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(A,2) - std(A,0,2)/sqrt(n),s)) ] ,...
    c(1,:), 'facealpha',0.3, 'edgecolor','none')
fill([erptime, fliplr(erptime)],...
    [smooth(mean(B,2) + std(B,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(B,2) - std(B,0,2)/sqrt(n),s)) ] ,...
    c(2,:), 'facealpha',0.3, 'edgecolor','none')


% [h,p ]= ttest(A',B');
% h(h ==0) = NaN;
% plot(erptime, h*1e-13, 'k*')
% ylim([-1 1]*1.5e-13); 
xlabel('time(s)'); ylabel('source dipole a.u.')
legend('Cue A','Cue B')

erfcon = [squeeze(ERF_AX(ix,:,:)); squeeze(ERF_AY(ix,:,:)) ; squeeze(ERF_BX(ix,:,:)); squeeze(ERF_BY(ix,:,:))];
vec = corr(erfcon);
% erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule
 subplot(122)
 
AX = squeeze(ERF_AX(ix,:,:)).*sign(vec(1,:));
AY = squeeze(ERF_AY(ix,:,:)).*sign(vec(1,:));
BX = squeeze(ERF_BX(ix,:,:)).*sign(vec(1,:));
BY = squeeze(ERF_BY(ix,:,:)).*sign(vec(1,:));

 
plot(erptime, smooth(mean(AX,2),s),'linewidth',2 )
hold all
plot(erptime,  smooth(mean(AY,2),s),'linewidth',2 )
plot(erptime,smooth( mean(BX,2),s) )
plot(erptime, smooth(mean(BY,2),s) )


fill([erptime, fliplr(erptime)],...
    [smooth(mean(AX,2) + std(AX,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(AX,2) - std(AX,0,2)/sqrt(n),s)) ] ,...
    c(1,:), 'facealpha',0.3, 'edgecolor','none')
fill([erptime, fliplr(erptime)],...
    [smooth(mean(AY,2) + std(AY,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(AY,2) - std(AY,0,2)/sqrt(n),s)) ] ,...
    c(2,:), 'facealpha',0.3, 'edgecolor','none')
% fill([erptime, fliplr(erptime)],...
%     [smooth(mean(BX,2) + std(BX,0,2)/sqrt(n),s) ; ...
%     flipud(smooth(mean(BX,2) - std(BX,0,2)/sqrt(n),s)) ] ,...
%     c(3,:), 'facealpha',0.3, 'edgecolor','none')
% fill([erptime, fliplr(erptime)],...
%     [smooth(mean(BY,2) + std(BY,0,2)/sqrt(n),s) ; ...
%     flipud(smooth(mean(BY,2) - std(BY,0,2)/sqrt(n),s)) ] ,...
%     c(4,:), 'facealpha',0.3, 'edgecolor','none')

% [h,p ]= ttest(AX',AY');
% h(h ==0) = NaN;
% plot(erptime, h*1e-13, 'k*')
grid on
% ylim([-1 1]*1.5e-13); 
xlabel('time(s)'); ylabel('source dipole a.u.')
legend('Probe AX','Probe AY','Probe BX','Probe BY')

% set(gcf,'name',sprintf('Corr ERF peak at %dms ',round(tt*1e3)))
% saveas(gcf,sprintf('%s/Comm_ERF_%dms_%s.jpg',resultsfolder,round(tt*1e3),stimname))

% title(sprintf('Incong Comm %d,%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm_%d-%dms_resp.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000))
% 
% saveas(gcf,sprintf('%s/Incong_Corr-base_%d-%dms_13-30Hz_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000,stimname))



%% Reactive control?

close all

s = 3;
% for tx = 0:0.05:0.4
%     iit = erptime >= tx & erptime <= (tx+0.1);
    iit = erptime >= 0.2 & erptime <= 0.35;

     % contrast with two conditions
    x1 = ( mean(ERF_AY(:,iit,s),2) );
    x2 = ( mean( ERF_AX(:,iit,s),2));

   
%     iitb = erptime >= -0.5 & erptime <= -0.1;
    mriplot = 'mni';
 
    T = zeros(sourcemodel.dim);
    T(sourcemodel.inside) = (abs(x1) - abs(x2)) * 1e15;
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



    crang = [-10 10];
    cfg = [];
    cfg.method        = 'ortho'; %'ortho'
    cfg.location   = 'max';

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

%     title(sprintf('B-A %d-%d ms',round(tx*1e3), round((tx+0.1)*1e3)))
    drawnow
% colormap default


%%
 iit = erptime >= 0.2 & erptime <= 0.5;

mnicor =  [-3 39 31]/10;
[xdiff,iix] = min(sqrt( sum((sourceant.pos - mnicor).^2 ,2) ) );

T = zeros(sourcemodel.dim);
T(iix) = 1;
T = T(sourcemodel.inside);
ix = find(T);

% erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule
figure; set(gcf,'color','w','position', [83   139   1575  433])
subplot(121)
A = squeeze(ERF_A(ix,:,s));
B = squeeze(ERF_B(ix,:,s));

AY = squeeze(ERF_AY(ix,:,s));
AX = squeeze(ERF_AX(ix,:,s));
plot(erptime, AX,'linewidth',2 )
hold on
plot(erptime, AY,'linewidth',2 )
xlim([0, 0.5])

hold all
plot(erptime,smooth( mean(B,2),s),'linewidth',2 )
grid on
c = get(gca,'colororder');
fill([erptime, fliplr(erptime)],...
    [smooth(mean(A,2) + std(A,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(A,2) - std(A,0,2)/sqrt(n),s)) ] ,...
    c(1,:), 'facealpha',0.3, 'edgecolor','none')
fill([erptime, fliplr(erptime)],...
    [smooth(mean(B,2) + std(B,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(B,2) - std(B,0,2)/sqrt(n),s)) ] ,...
    c(2,:), 'facealpha',0.3, 'edgecolor','none')


% [h,p ]= ttest(A',B');
% h(h ==0) = NaN;
% plot(erptime, h*1e-13, 'k*')
% ylim([-1 1]*1.5e-13); 
xlabel('time(s)'); ylabel('source dipole a.u.')
legend('Cue A','Cue B')

erfcon = [squeeze(ERF_AX(ix,:,:)); squeeze(ERF_AY(ix,:,:)) ; squeeze(ERF_BX(ix,:,:)); squeeze(ERF_BY(ix,:,:))];
vec = corr(erfcon);
% erfcon = erfcon.*sign(vec(1,:)); % aligns to the left precentral lobule
 subplot(122)
 
AX = squeeze(ERF_AX(ix,:,:)).*sign(vec(1,:));
AY = squeeze(ERF_AY(ix,:,:)).*sign(vec(1,:));
BX = squeeze(ERF_BX(ix,:,:)).*sign(vec(1,:));
BY = squeeze(ERF_BY(ix,:,:)).*sign(vec(1,:));

 
plot(erptime, smooth(mean(AX,2),s),'linewidth',2 )
hold all
plot(erptime,  smooth(mean(AY,2),s),'linewidth',2 )
plot(erptime,smooth( mean(BX,2),s) )
plot(erptime, smooth(mean(BY,2),s) )


fill([erptime, fliplr(erptime)],...
    [smooth(mean(AX,2) + std(AX,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(AX,2) - std(AX,0,2)/sqrt(n),s)) ] ,...
    c(1,:), 'facealpha',0.3, 'edgecolor','none')
fill([erptime, fliplr(erptime)],...
    [smooth(mean(AY,2) + std(AY,0,2)/sqrt(n),s) ; ...
    flipud(smooth(mean(AY,2) - std(AY,0,2)/sqrt(n),s)) ] ,...
    c(2,:), 'facealpha',0.3, 'edgecolor','none')
% fill([erptime, fliplr(erptime)],...
%     [smooth(mean(BX,2) + std(BX,0,2)/sqrt(n),s) ; ...
%     flipud(smooth(mean(BX,2) - std(BX,0,2)/sqrt(n),s)) ] ,...
%     c(3,:), 'facealpha',0.3, 'edgecolor','none')
% fill([erptime, fliplr(erptime)],...
%     [smooth(mean(BY,2) + std(BY,0,2)/sqrt(n),s) ; ...
%     flipud(smooth(mean(BY,2) - std(BY,0,2)/sqrt(n),s)) ] ,...
%     c(4,:), 'facealpha',0.3, 'edgecolor','none')

% [h,p ]= ttest(AX',AY');
% h(h ==0) = NaN;
% plot(erptime, h*1e-13, 'k*')
grid on
% ylim([-1 1]*1.5e-13); 
xlabel('time(s)'); ylabel('source dipole a.u.')
legend('Probe AX','Probe AY','Probe BX','Probe BY')

% set(gcf,'name',sprintf('Corr ERF peak at %dms ',round(tt*1e3)))
% saveas(gcf,sprintf('%s/Comm_ERF_%dms_%s.jpg',resultsfolder,round(tt*1e3),stimname))

% title(sprintf('Incong Comm %d,%dms',timew(1)*1e3,timew(2)*1e3))
% saveas(gcf,sprintf('%s/Incong_absComm_%d-%dms_resp.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000))
% 
% saveas(gcf,sprintf('%s/Incong_Corr-base_%d-%dms_13-30Hz_%s.jpg',resultsfolder,...
%     timew(1)*1000, timew(2)*1000,stimname))


%% Plot Freq

crang = [-4 4];

mriplot = 'mni';
plotopt = 'slice';  % 'slice'; %'ortho'

T = zeros(sourcemodel.dim);
T(sourcemodel.inside) = mean(Theta_B,2) -  mean(Theta_A,2) ; 
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