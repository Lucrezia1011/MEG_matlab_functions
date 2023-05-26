function [buttonsamp,rt] = match_responses(stimsamp, buttonresp, resp, f) 

buttonsamp = zeros(size(stimsamp ));
rt  = zeros(size(stimsamp ));
for tt = 1:length(stimsamp ) % left button press to probe
    bsamp = buttonresp.(resp)( (buttonresp.(resp) - stimsamp(tt))>0 & (buttonresp.(resp) - stimsamp(tt))< f);
    if ~isempty(bsamp)  % only use first button press
        buttonsamp(tt) = bsamp(1);
        rt(tt) = (buttonsamp(tt) - stimsamp(tt))/f;
    end
    
end
% buttonsamp(buttonsamp == 0) = [];