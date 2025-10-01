function [cuesamp,respsamp,rt,resp_bool]= flanker_trials(cuesamp,buttonpress,resp,f,minrt,maxrt)
  
    [respsamp,rt_corr] = match_responses(cuesamp, buttonpress,resp , f,minrt,maxrt);
    if strcmp(resp,'right')
        respcomm = 'left';
    else
        respcomm = 'right';
    end
    [respsamp_comm,rt_comm] = match_responses(cuesamp, buttonpress, respcomm, f,minrt,maxrt);
    resp_bool = rt_corr > 0;
    resp_bool = resp_bool(rt_corr | rt_comm); % exclude omission trials
    cuesamp = cuesamp(respsamp | respsamp_comm);  % exclude omission trials
    respsamp(rt_comm>0) = respsamp_comm(rt_comm>0); 
    respsamp = respsamp(respsamp | respsamp_comm);  % exclude omission trials
    rt = rt_corr;
    rt(rt_comm>0) = rt_comm(rt_comm>0);
    rt = rt(rt_corr | rt_comm); % exclude omission trials

   