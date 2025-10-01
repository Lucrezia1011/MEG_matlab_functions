function [buttonsamp,flansamp] = resp2trial(resp,cuesamp,buttonresp,f)
if strcmp(resp,'left')
    resp2 = 'right';
else
    resp2 = 'left';
end
buttonsamp = zeros(size(cuesamp));
flansamp= zeros(size(cuesamp));
for tti = 1:length(cuesamp) % left button press to cue
    bsamp = buttonresp.(resp)( (buttonresp.(resp) - cuesamp(tti))>0 & ...
        (buttonresp.(resp) - cuesamp(tti))< f & ...
        ~any(((buttonresp.(resp) - buttonresp.(resp2)')>0 & (buttonresp.(resp) - buttonresp.(resp2)')<f),2));
    if ~isempty(bsamp)  % only use first button press
        buttonsamp(tti) = bsamp(1);
        flansamp(tti) = cuesamp(tti);
    end
end
buttonsamp(buttonsamp == 0) = [];
flansamp(flansamp == 0) = [];
end