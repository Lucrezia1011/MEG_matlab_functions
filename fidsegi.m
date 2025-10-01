function elec_new = fidsegi(mri,fids_name,elecfilename,plotOpt)
% Lucrezia Liuzzi, last updated 2023/12/07
% Transform nifti mri into ctf coordinates
% mriC = fids2ctf(mri_name,fids_name,plotOpt)
%
% mri_name :  file name of MRI in nift format
% fids_name : fiducial position .tag file, can be left blank if fiducial
% coordinates are written in mri .json file
% plotOpt  : =1 for plotting

fids_inds = zeros(3,4);
fids_inds(:,4) = 1;

fileID = fopen(fids_name,'r');
fids_char = fscanf(fileID,'%c');
fclose(fileID);

% 4 lines of 66 characters each
for iF = 1:3 % transform to RAS
    fids_inds(iF,1) = -str2double(fids_char(66*iF+(18:28))); 
    fids_inds(iF,2) = -str2double(fids_char(66*iF+(30:40)));
    fids_inds(iF,3) = str2double(fids_char(66*iF+(42:52)));
end

% transform afni fiducial tags to voxel coords
fid_vox = round(inv(mri.transformorig)*fids_inds')';
% transform voxel coords to current ctf space
fids_ctf = fid_vox*mri.transform'; % 

mri.elecpos = [
    fids_ctf(1,1:3); %position of nasion
    fids_ctf(2,1:3); %position of LPA
    fids_ctf(3,1:3); %position of RPA
    ];
mri.label = {'Nasion', 'LPA', 'RPA'};


% Check with plots ---------------------------------------------------
if plotOpt == 1
    fids_labels = {'Nasion';'Left Ear';'Right Ear'};
    figure;
    clf
    for iF = 1:3
        subplot(3,3,(iF-1)*3+1)
        imagesc(rot90(squeeze(mri.anatomy(fid_vox(iF,1),:,:)))); axis equal;
        set(gca,'Xdir','reverse')
        title(fids_labels{iF})
        subplot(3,3,(iF-1)*3+2)
        imagesc(rot90(squeeze(mri.anatomy(:,fid_vox(iF,2),:)))); axis equal
        subplot(3,3,(iF-1)*3+3)
        imagesc(fliplr(mri.anatomy(:,:,fid_vox(iF,3)))'); axis equal

    end
    colormap gray
end
% --------------------------------------------------------------------
geo = readtable(elecfilename);
geo(:,5) = [];
geo = renamevars(geo,["Var1","Var2", "Var3", "Var4"],["labels","x","y","z"]);
%Var2: post negative - front positive
%Var3: left negative - right positive
%Var4: infer negative - superior positive
elec = [];
elec.pos = [geo.x, -geo.y, geo.z];
elec.label = geo.labels;
% elec.coordsys = 'ctf';
elec.unit = 'cm';
elec = ft_convert_units(elec,'mm'); % should be the same unit as MRI

elec.fid = [];
elec.fid.unit = 'mm'; 
elec.fid.elecpos = [geo.x(1:3), -geo.y(1:3), geo.z(1:3)]*10;
elec.fid.label = {'Nasion';'RPA';'LPA'  };
% ----------------------------------------------------------
% coregister the electrodes to the MRI using fiducials

cfg = [];
cfg.fiducial = {'Nasion';'RPA';'LPA'  };
cfg.method = 'fiducial';
cfg.target = mri;
elec_new = ft_electroderealign(cfg,elec);


% cfg = [];
% % cfg.fiducial = {'Nasion';'RPA';'LPA'  };
% cfg.method = 'fiducial';
% cfg.fiducial.nas    = fid_vox(1,1:3); %position of nasion
% cfg.fiducial.lpa    = fid_vox(2,1:3); %position of LPA
% cfg.fiducial.rpa    = fid_vox(3,1:3); %position of RPA
% cfg.viewresult     = 'yes';
% cfg.coordsys = 'mni';
% mriC = ft_volumerealign(cfg,mri);

% 
% figure; plot3(elec.pos(:,1),elec.pos(:,2),elec.pos(:,3),'k.')
% axis equal; hold on
% plot3(elec.fid.elecpos(1,1), elec.fid.elecpos(1,2), elec.fid.elecpos(1,3),'b.','MarkerSize',15); 
% plot3(elec.fid.elecpos(2,1), elec.fid.elecpos(2,2), elec.fid.elecpos(2,3),'r.','MarkerSize',15); 
% plot3(elec.fid.elecpos(3,1), elec.fid.elecpos(3,2), elec.fid.elecpos(3,3),'g.','MarkerSize',15); 
% 
% 
% figure; plot3(elec_new.chanpos(:,1),elec_new.chanpos(:,2),elec_new.chanpos(:,3),'k.')
% axis equal


