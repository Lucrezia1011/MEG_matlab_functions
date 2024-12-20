function grid = mniLeadfields_multiSpheres(data_name,processing_folder,gridres,mri, varargin)
% Calculate lead fields with fieldtrip with local spheres approximation on 
% regularly spaced MNI grid warped to individual anatomy 
% 
% grid = mniLeadfields_multiSpheres(data_name,processing_folder,gridres,mri)
% data_name          = dataset name (.ds)
% processing_folder  = folder for data derivatives
% gridres            = beamformer grid resolution in mm
% mri                = co-registered mri
% keepbad            = keep bad channels (for MEG realigned or other channels interpolation)

if isempty(varargin)
    keepbad = 0;
else
    keepbad = varargin{1};
end

leadfield_name =sprintf( '%s%s/leadfields_multiSpheres_%.0fmm.mat',processing_folder,data_name(1:end-3),gridres);
if ~exist(leadfield_name,'file')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Segment MRI
    sens = ft_read_sens(data_name,'senstype','meg');


    if ~exist([processing_folder,'/headmodel_multiSpheres.mat'],'file')
        cfg = [];
        cfg.output  = 'brain';
        segmentmri = ft_volumesegment(cfg,mri);
        
        segmentmri.anatomy = mri.anatomy;
        
        % Plot mri and brain volume for debugging
        cfg = [];
        cfg.anaparameter = 'anatomy';
        cfg.funparameter = 'brain';
        cfg.location = [0 0 60];
        ft_sourceplot(cfg, segmentmri)
        
        %         % no ref sensors
        %         grad = sens;
        %         sensgrad= strcmp(sens.chantype,'meggrad');
        %         grad.chanori(~sensgrad,:) = [];
        %         grad.chanpos(~sensgrad,:) = [];
        %         grad.chantype(~sensgrad) = [];
        %         grad.chanunit(~sensgrad) = [];
        %         grad.label(~sensgrad) = [];
        %         grad.tra(~sensgrad,:) = [];
        
%         cfg             = [];
%         cfg.tissue      = {'brain'};
%         cfg.numvertices = [2400];
%         bnd    = ft_prepare_mesh(cfg,segmentmri);
        
        cfg = [];
        cfg.method = 'localspheres';%'singleshell';
        cfg.grad            = sens;
        cfg.numvertices =  6000; % increase number of verticies from default 3000
        %         cfg.radius  = 150;  % 85mm default, 4 spikes at top, using 70mm made it worse
        cfg.feedback = 'no';
        vol = ft_prepare_headmodel(cfg, segmentmri);
%         
%         figure
%         ft_plot_mesh(bnd,'unit','cm','facealpha',.5); hold on
%         ft_plot_sens(sens, 'unit', 'cm'); hold on;
%         ft_plot_headmodel(vol, 'facecolor', 'cortex','edgecolor',[0,0,0], 'grad', sens, 'unit', 'cm','facealpha','0.5');
%         title(['radius ',num2str(cfg.radius),'mm'])
       
       
%         figure; cla
%         ft_plot_mesh(bnd,'unit','mm','facealpha',.5); hold on
% 
%        view([-20 20])
%          headmodel = ft_convert_units(vol, sens.unit);
%         [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens);
%         [bnd.pos, bnd.tri] = headsurface(headmodel, sens,'surface','brain');
        
        save([processing_folder,'/headmodel_multiSpheres.mat'],'vol')
    else
        load([processing_folder,'/headmodel_multiSpheres.mat']);
        if length(vol.label) < length(sens.label)
            cfg = [];
            cfg.output  = 'brain';
            segmentmri = ft_volumesegment(cfg,mri);
            
            segmentmri.anatomy = mri.anatomy;
            
            % Plot mri and brain volume for debugging
            cfg = [];
            cfg.anaparameter = 'anatomy';
            cfg.funparameter = 'brain';
            cfg.location = [0 0 60];
            ft_sourceplot(cfg, segmentmri)

            cfg = [];
            cfg.method = 'localspheres';%'singleshell';
            cfg.grad            = sens;
            cfg.numvertices =  6000; % increase number of verticies from default 3000
            %         cfg.radius  = 150;  % 85mm default, 4 spikes at top, using 70mm made it worse
            cfg.feedback = 'no';
            vol = ft_prepare_headmodel(cfg, segmentmri);
            save([processing_folder,'/headmodel_multiSpheres.mat'],'vol')
        end
    
    end
    
    
   
    %% MNI template brain
    
    % Load fieldtrip MNI grid
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    template_grid = sourcemodel;
    
    clear sourcemodel
    
    %% Sourcemodel warp MNI grid
    
    % sourcemodel based on 5mm grid MNI brain
    cfg = [];
    cfg.mri = mri;
    cfg.warpmni = 'yes';
    cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
    cfg.unit      = 'm';
    cfg.nonlinear = 'yes';
    sourcemodel = ft_prepare_sourcemodel(cfg);
    locs = sourcemodel.pos;
   
    %% Calculate lead fields
      
    cfg                 = [];
    cfg.grad            = sens;
    cfg.headmodel       = vol;
    cfg.reducerank      = 2;
    cfg.channel         = {'MEG'};
%     cfg.sourcemodel.pos = locs; %sourcemodel.pos
%     cfg.sourcemodel.unit   = 'm';
    cfg.sourcemodel    = sourcemodel; % Thia works
    cfg.siunits         = true;
    cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
    [grid] = ft_prepare_leadfield(cfg);
    
    %% Eliminate Bad channels    
    if keepbad == 0
        % Get Bad channel names
        fid = fopen([data_name,'/BadChannels']);
        BadChannels = textscan(fid,'%s');
        fclose(fid);
    
        % Delete Bad channels
        chanInd = zeros(size(grid.label));
        for iiC = 1:length(BadChannels{1})
            chanInd = chanInd | strcmp(grid.label,BadChannels{1}{iiC});
        end
        grid.label(find(chanInd)) = [];
        for ii = find(grid.inside)'
            grid.leadfield{ii}((find(chanInd)),:) = [];
        end
    end
    %%
    
    save(leadfield_name,'grid');
else
    load(leadfield_name);
end
