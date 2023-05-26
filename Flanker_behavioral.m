subs = {'24531';'24580';'24590';'24584'};

nerrs = zeros(1,length(subs));
nmisses = zeros(1,length(subs));
for ss = 1:length(subs)
    sub = subs{ss};
    datapath = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/meg/'];

    if ss == 1
        datafolder = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/ses-02/behavioral/'];
    else
        datafolder = ['/data/liuzzil2/MEG_AXCPT_Flanker/data/sub-',sub,'/behavioral/'];
    end
    
    d = dir(datafolder);
    
    bnames = cell(1,2);
    n = 0;
    nerr = 0;
    nmiss = 0;
    for ii = 3:length(d)
        if contains(d(ii).name , 'Flanker3') && contains(d(ii).name , '.csv')
            n  = n+1;
            bnames{n} =  d(ii).name;
            bv = readtable([datafolder,bnames{n}]);
            nerr = nerr + nnz( strcmp(bv.key_resp_resp_type ,'commission') ) ;
            nmiss = nmiss +  nnz( strcmp(bv.key_resp_resp_type ,'omission') ) ;
        end
    end
   
    nerrs(ss) = nerr;
    nmisses(ss) = nmiss;
    % 192 * 2 trials
    
end