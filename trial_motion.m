function cc_dem = trial_motion(headpos,trials)
% trials = vector with trial numbers, e.g. 1:ntrials

ntrials = length(trials);
% calculate the mean coil position per trial
coil1 = zeros(3,ntrials);
coil2 = zeros(3,ntrials);
coil3 = zeros(3,ntrials);
for ii = 1:ntrials
    t = trials(ii);
    coil1(:,ii) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
    coil2(:,ii) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
    coil3(:,ii) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
end

% calculate the headposition and orientation per trial
cc = circumcenter(coil1, coil2, coil3);

% demean to obtain translations and rotations from the average position and orientation
% transpose to construct a nsamples-by-nregressors matrix
cc_dem = [cc - repmat(mean(cc,2),1,size(cc,2))]';