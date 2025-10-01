function [buttonsamp,rt] = match_responses(stimsamp, buttonresp, resp, f, varargin ) 
if nargin == 6
    minrt = varargin{1};
    maxrt = varargin{2};
    minrt = minrt*f; maxrt = maxrt*f;
elseif nargin == 4 % old definition  
    minrt = 0;
    maxrt = f; % 1s
end

if strcmp(resp,'left')
    resp2 = 'right';
else
    resp2 = 'left';
end

buttonsamp = zeros(size(stimsamp ));
rt  = zeros(size(stimsamp ));
for tt = 1:length(stimsamp ) % left button press to probe
    bsamp = buttonresp.(resp)( (buttonresp.(resp) - stimsamp(tt))>=minrt & (buttonresp.(resp) - stimsamp(tt))<= maxrt);
    bsamp2 = buttonresp.(resp2)( (buttonresp.(resp2) - stimsamp(tt))>=minrt & (buttonresp.(resp2) - stimsamp(tt))<= maxrt);
    if ~isempty(bsamp)  % only use first button press
        if isempty(bsamp2) % no other button was pressed
            buttonsamp(tt) = bsamp(1);
            rt(tt) = (buttonsamp(tt) - stimsamp(tt))/f;
        else % other button was pressed
            if bsamp2(1) > bsamp(1) % incorrect button was pressed after correct 
                buttonsamp(tt) = bsamp(1);
                rt(tt) = (buttonsamp(tt) - stimsamp(tt))/f;
            end
        end
    end
    
end
% buttonsamp(buttonsamp == 0) = [];