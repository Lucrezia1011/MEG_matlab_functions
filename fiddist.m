function fiddist(sdan, ses, datapath)
% fiddist(sdan, ses, datapath)
% Lucrezia Liuzzi, last updated 2024/01/23
% equivalent of fiddist.py

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


datapathsub = [datapath,'sub-',sdan,'/ses-',ses,'/'];
fids_name = [datapathsub,'meg/sub-',sdan,'_fiducials.tag'];
fids_inds = zeros(3,3);


if ~exist(fids_name,'file')
    fids_name = [datapathsub,'meg/sub-',sdan,'_ses-',ses,'_fiducials.tag'];
end
if ~exist(fids_name,'file')
    error('No fiducials file!')
end

fileID = fopen(fids_name,'r');
fids_char = fscanf(fileID,'%c');
fclose(fileID);

% 4 lines of 66 characters each
for iF = 1:3
    fids_inds(iF,1) = str2double(fids_char(66*iF+(18:28))); 
    fids_inds(iF,2) = str2double(fids_char(66*iF+(30:40)));
    fids_inds(iF,3) = str2double(fids_char(66*iF+(42:52)));
end


% % Read fiducial coordinates from  mri json file
% json_name = [datapathsub,'meg/sub-',sdan,'_coordsystem.json'];
% 
% ft_info('reading %s\n', json_name);
% ft_hastoolbox('jsonlab', 1);
% json = loadjson(json_name);

megdir = dir([datapathsub,'meg/']);
nfiles = length(megdir)-5;
datanames = cell(nfiles,1);
fids_meg = zeros(3,3,nfiles);
n = 0;
for jj = 3:length(megdir)
    if strcmp(megdir(jj).name(end-2:end), '.ds')
        n = n+1;
        datanames{n} = megdir(jj).name;
        hdr = ft_read_header([datapathsub,'meg/',datanames{n}]);
        fids_meg(:,:,n) = hdr.orig.hc.head'; 
        
    elseif contains(megdir(jj).name, 'T1w+orig')
        mriname = megdir(jj).name;
        mrinamer = mriname(1:strfind(mriname,'+')-1);
    end
end
fids_meg = fids_meg(:,:,1:n);
fids_meg_av = mean(fids_meg,3);

% clc
% % disp('MEG fiducial coordinates:')
% % disp(fids_meg)
% disp('Average MEG fiducial coordinates:')
% disp(fids_meg_av)
% disp('stdev of MEG fiducial coordinates:')
% disp(std(fids_meg,0,3))
% 
% 
% %%
% 
% 
% dlr_mri = sqrt(sum((fids_inds(2,:) - fids_inds(3,:)).^2));
% dlr_meg = sqrt(sum((fids_meg_av(2,:) - fids_meg_av(3,:)).^2))*10;
% 
% fprintf('Total Left and right ear distance:\nMEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
%     dlr_meg,dlr_mri,(dlr_mri-dlr_meg),abs(dlr_mri-dlr_meg)<5)
% 
% dnl_mri = sqrt(sum((fids_inds(1,:) - fids_inds(2,:)).^2));
% dnl_meg = sqrt(sum((fids_meg_av(1,:) - fids_meg_av(2,:)).^2))*10;
% 
% fprintf('Total nasion and left ear distance:\nMEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
%     dnl_meg,dnl_mri,(dnl_mri-dnl_meg),abs(dnl_mri-dnl_meg)<5)
% 
% dnr_mri = sqrt(sum((fids_inds(1,:) - fids_inds(3,:)).^2));
% dnr_meg = sqrt(sum((fids_meg_av(1,:) - fids_meg_av(3,:)).^2))*10;
% 
% fprintf('Total nasion and right ear distance:\nMEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
%     dnr_meg,dnr_mri,(dnr_mri-dnr_meg),abs(dnr_mri-dnr_meg)<5)

%% transformation matrix


mri_name = [ datapathsub, '/anat/',mrinamer,'.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
mri = ft_read_mri(mri_name,'dataformat','nifti');

T = mri.transform; % = mri.hdr.vox2ras1, shifted by one voxel 
if sign(T(1,1)) == 1 && sign(T(2,2)) == 1 && sign(T(3,3)) == 1
    T(1:2,4) = -T(1:2,4); % adjust transformation matrix from fsl to afni (LPS)
    T(1,1) = -1; T(2,2) = -1;
elseif  sign(T(1,3)) == -1 && sign(T(2,1)) == -1 && sign(T(3,2)) == -1
    T(1,3) = -T(1,3); T(2,1) = -T(2,1);
    T(1:2,4) = -T(1:2,4); % adjust transformation matrix from RAS to afni (LPS)

else
    error('Unrecognized MRI transform')
end

fids_inds_flip =  zeros(3,4);
fids_inds_flip(:,1:3) = fids_inds;
fids_inds_flip(:,4) = 1;
% MRI voxels
fid_vox = round(inv(T)*fids_inds_flip')';


cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fid_vox(1,1:3); %position of nasion
cfg.fiducial.lpa    = fid_vox(2,1:3); %position of LPA
cfg.fiducial.rpa    = fid_vox(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

mriC = ft_volumerealign(cfg,mri);
mri_vox = mriC.transform*fid_vox';
mri_vox = mri_vox(1:3,:)';

%%

clc
% disp('MEG fiducial coordinates:')
% disp(fids_meg)
disp('Average MEG fiducial coordinates:')
disp('    PA         RL         IS')
disp(fids_meg_av*10)
disp('MRI fiducial coordinates:')
disp(mri_vox)

disp('MEG - MRI:')
disp(fids_meg_av*10 - mri_vox)



dlr_mri = sqrt(sum((mri_vox(2,:) - mri_vox(3,:)).^2));
dlr_meg = sqrt(sum((fids_meg_av(2,:) - fids_meg_av(3,:)).^2))*10;

fprintf('Total Left and right ear distance:\nMEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
    dlr_meg,dlr_mri,(dlr_mri-dlr_meg),abs(dlr_mri-dlr_meg)<5)

dnl_mri = sqrt(sum((mri_vox(1,:) - mri_vox(2,:)).^2));
dnl_meg = sqrt(sum((fids_meg_av(1,:) - fids_meg_av(2,:)).^2))*10;

fprintf('Total nasion and left ear distance:\nMEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
    dnl_meg,dnl_mri,(dnl_mri-dnl_meg),abs(dnl_mri-dnl_meg)<5)

dnr_mri = sqrt(sum((mri_vox(1,:) - mri_vox(3,:)).^2));
dnr_meg = sqrt(sum((fids_meg_av(1,:) - fids_meg_av(3,:)).^2))*10;

fprintf('Total nasion and right ear distance:\nMEG %.1fmm\nMRI %.1fmm\ndiff = %.1fmm good=%d\n\n',...
    dnr_meg,dnr_mri,(dnr_mri-dnr_meg),abs(dnr_mri-dnr_meg)<5)

%%
% dx_lr_mri = mri_vox(2,2) - mri_vox(3,2);
% dx_lr_meg = (fids_meg_av(2,2) - fids_meg_av(3,2))*10;
% fprintf('dx Left and right ear distance:\nMRI %.1fmm\nMEG %.1fmm\n\n',dx_lr_mri,dx_lr_meg)
% 
% dyz_mri = mri_vox(1,1);
% dyz_meg = fids_meg_av(1,1)*10;
% fprintf('dyz nasion to mean ear distance:\nMRI %.1fmm\nMEG %.1fmm\n\n',dyz_mri,dyz_meg)


%%
% clc
% dx_lr_mri = fids_inds(2,1) - fids_inds(3,1);
% dx_lr_meg = (json.HeadCoilCoordinates.left_ear(2) - json.HeadCoilCoordinates.right_ear(2))*10;
% fprintf('dx Left and right ear distance:\nMRI %dmm\nMEG %.0fmm\n\n',dx_lr_mri,dx_lr_meg)
% 
% dlr_mri = sqrt(sum((fids_inds(2,:) - fids_inds(3,:)).^2));
% dlr_meg = sqrt(sum((json.HeadCoilCoordinates.left_ear - json.HeadCoilCoordinates.right_ear).^2))*10;
% 
% fprintf('Total Left and right ear distance:\nMRI %.1fmm\nMEG %.1fmm\ndiff = %.1fmm good=%d\n\n',...
%     dlr_mri,dlr_meg,abs(dlr_mri-dlr_meg),abs(dlr_mri-dlr_meg)<5)
% 
% dnl_mri = sqrt(sum((fids_inds(1,:) - fids_inds(2,:)).^2));
% dnl_meg = sqrt(sum((json.HeadCoilCoordinates.nasion - json.HeadCoilCoordinates.left_ear).^2))*10;
% 
% fprintf('Total nasion and left ear distance:\nMRI %.1fmm\nMEG %.1fmm\ndiff = %.1fmm good=%d\n\n',...
%     dnl_mri,dnl_meg,abs(dnl_mri-dnl_meg),abs(dnl_mri-dnl_meg)<5)
% 
% 
% dnr_mri = sqrt(sum((fids_inds(1,:) - fids_inds(3,:)).^2));
% dnr_meg = sqrt(sum((json.HeadCoilCoordinates.nasion - json.HeadCoilCoordinates.right_ear).^2))*10;
% 
% fprintf('Total nasion and right ear distance:\nMRI %.1fmm\nMEG %.1fmm\ndiff = %.1fmm good=%d\n\n',...
%     dnr_mri,dnr_meg,abs(dnr_mri-dnr_meg),abs(dnr_mri-dnr_meg)<5)
% 
% dyz_mri = sqrt((fids_inds(1,3) - mean(fids_inds(2:3,3)))^2 + ...
%     (fids_inds(1,2) - mean(fids_inds(2:3,2)))^2);
% dyz_meg = json.HeadCoilCoordinates.nasion(1)*10;
% fprintf('dyz nasion to mean ear distance:\nMRI %.0fmm\nMEG %.0fmm\n\n',dyz_mri,dyz_meg)
