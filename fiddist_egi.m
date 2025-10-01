function fiddist_egi(sdan, ses, datapath)
% fiddist(sdan, ses, datapath)
% Lucrezia Liuzzi, last updated 2024/01/23
% equivalent of fiddist.py
if ~exist('fieldtrip2ctf.m','file')
    addpath /data/liuzzil2/fieldtrip/
    ft_defaults
end
datapathsub = [datapath,'/BIDS/sub-',sdan,'/ses-',ses,'/'];
eegdir = [datapathsub,'eeg/'];
derpathsub = [datapath,'/derivatives/sub-',sdan,'/'];
fids_name = [derpathsub,'/sub-',sdan,'_ses-',ses,'_egimarkers.tag'];
d = dir(eegdir);
elecfilename = d(startsWith({d.name},'GeoScan')).name;
d = dir(derpathsub);
mriname = d(endsWith({d.name},'T1w+orig.BRIK')).name;
mrinamer = mriname(1:strfind(mriname,'+')-1);
 
mri_name = [ datapathsub, '/anat/',mrinamer,'.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
mri = ft_read_mri(mri_name,'dataformat','nifti');


if ~exist(fids_name,'file')
    error('No fiducials file!')
end

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
fid_vox = round(inv(mri.transform)*fids_inds')';



% --------------------------------------------------------------------
geo = readtable([eegdir,elecfilename]);
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
elec.fid.elecpos = [geo.x(contains(elec.label,'Fid')), -geo.y(contains(elec.label,'Fid')), geo.z(contains(elec.label,'Fid'))]*10;
elec.fid.label = {'Nasion';'RPA';'LPA'  };
% --------------------------------------------------------------------

cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fid_vox(1,1:3); %position of nasion
cfg.fiducial.lpa    = fid_vox(2,1:3); %position of LPA
cfg.fiducial.rpa    = fid_vox(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

mriC = ft_volumerealign(cfg,mri);

fids_ctf = fid_vox*mriC.transform'; % 

mriC.elecpos = [
    fids_ctf(1,1:3); %position of nasion
    fids_ctf(2,1:3); %position of LPA
    fids_ctf(3,1:3); %position of RPA
    ];
mriC.label = {'Nasion', 'LPA', 'RPA'};

cfg = [];
cfg.fiducial = {'Nasion';'RPA';'LPA'  };
cfg.method = 'fiducial';
cfg.target = mriC;
elec_new = ft_electroderealign(cfg,elec);


%%
eeg_elec = elec_new.fid.elecpos([1,3,2],:);
mri_vox = fids_ctf(:,1:3);
clc
% disp('MEG fiducial coordinates:')
% disp(fids_meg)
disp('Average EEG fiducial coordinates:')
disp('    PA         RL         IS')
disp(elec_new.fid.elecpos)
disp('MRI fiducial coordinates:')
disp(mri_vox)

disp('EEG - MRI:')
disp(eeg_elec - mri_vox)



dlr_mri = sqrt(sum((mri_vox(2,:) - mri_vox(3,:)).^2));
dlr_meg = sqrt(sum((eeg_elec(2,:) - eeg_elec(3,:)).^2));

fprintf('Total Left and right ear distance:\nEEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
    dlr_meg,dlr_mri,(dlr_mri-dlr_meg),abs(dlr_mri-dlr_meg)<5)

dnl_mri = sqrt(sum((mri_vox(1,:) - mri_vox(2,:)).^2));
dnl_meg = sqrt(sum((eeg_elec(1,:) - eeg_elec(2,:)).^2));

fprintf('Total nasion and left ear distance:\nEEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
    dnl_meg,dnl_mri,(dnl_mri-dnl_meg),abs(dnl_mri-dnl_meg)<5)

dnr_mri = sqrt(sum((mri_vox(1,:) - mri_vox(3,:)).^2));
dnr_meg = sqrt(sum((eeg_elec(1,:) - eeg_elec(3,:)).^2));

fprintf('Total nasion and right ear distance:\nEEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
    dnr_meg,dnr_mri,(dnr_mri-dnr_meg),abs(dnr_mri-dnr_meg)<5)

