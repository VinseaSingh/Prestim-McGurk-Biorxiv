function [GlobCoh] = Function_TmAvgGlobCoh(data, win, params, cfg)

%% Function to compute the Time-Averaged Global Coherence for a group (based on the percentage of illusory percept) of participants

%{
 Inputs:
    data - The preprocessed data that contains subjectwise EEG data 
           (a 1D cell array, each cell comprising data from a single subject)
    win     - Duration of the non-overlapping window (for the use in the Chronux function CrossSpecMatc.m)
    params  - structure with fields tapers, pad, Fs, fpass (for the use in the Chronux function CrossSpecMatc.m)
    cfg     -
           cfg.choice   - 'InterSub' for group-wise(based on percentage of /ta/) and subjectwise global coherence
                          'InterTrial_Groupwise' for global coherence based on the responses (/ta/ and /pa/) grouped over all participants
           cfg.interval - the size of the interval(%age of /ta/) for
                          grouping the participants e.g.- 20, 25, 50 etc.

 Outputs:
    GlobCoh  - A structured file with the outputs of function CrossSpecMatc.m
               for trails sorted based on percentage of /ta/ responses or based on type of responses (/ta/ and /pa/)
               Sc (cross spectral matrix frequency x channels x channels)
               Cmat Coherence matrix frequency x channels x channels
               Ctot Total coherence: SV(1)^2/sum(SV^2) (frequency)
               f (frequencies)

    GlobCoh.Trialnum - The trial indx from which the GlobCoh has been computed
                     - for choice 'Intersub'- Trialnum is a strc array that consists
                       of a cell array field 'participant-wise' comprising the indices of the trials used from each participant and
                       field 'Groupwise' that comprise the indices of the trials used from each group.

                     - for choice 'InterTrial'- is an array comprising the indices of the sorted trials (/ta/ or /pa/) higher in number that
                       that has been used in the computation of global coherence
%}


% Grouping participants based on their %age of /ta/ percept and computing the respective Global Coherence

choice = cfg.choice;

switch choice
    
    % In this case we pick equal number of trials from each participant that is
    % equivalent to the minimum trial numbers among all participants.
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        
    case 'InterSub'
        
        % extracting the behavior(/ta/ and /pa/ %ages) and trials from each subject
        
        Sub_len = length(data);            % number of subjects
        Percent_McG_ta = zeros(Sub_len,1);    % stores the %age of /ta/ response of each subject
        Percent_McG_pa = zeros(Sub_len,1);    % stores the %age of /pa/ response of each subject
        trnum_ta  = zeros(Sub_len,1);         % stores the NUMBER of /ta/ trials of each subject
        trnum_pa  = zeros(Sub_len,1);         % stores the NUMBER of /pa/ trials of each subject
        trnum_tot = zeros(Sub_len,1);         % stores the total NUMBER of trials (/ta/+/pa/)  of each subject
        
        for sub = 1:length(data)
            Percent_McG_ta(sub) = data{sub}.percent_McG.ta;
            Percent_McG_pa(sub) = data{sub}.percent_McG.pa;
            
            trnum_ta(sub) = size(data{sub}.McG_ta,3);
            trnum_pa(sub) = size(data{sub}.McG_pa,3);                          
        end
        trnum_tot = trnum_ta + trnum_pa;    % array that stores the total number of trials(ta+pa) of each subject
        
        interval = cfg.interval;
        Grpnum = 100/interval;              % number of groups
        Subject_Group = cell(Grpnum,1);     % stores the index of the subjects
        minPrcnt = 0; maxPrcnt = interval;
        
        for grp = 1: Grpnum          
            Subject_Group{grp}.Subindx = find(Percent_McG_ta >= minPrcnt & Percent_McG_ta < maxPrcnt );           
            if grp ~= length(Subject_Group) - 1
                minPrcnt = minPrcnt + interval;
                maxPrcnt = maxPrcnt + interval;
            else
                minPrcnt = minPrcnt + interval;
                maxPrcnt = maxPrcnt + interval + 1;  % this makes sure we include the participants with 100% ta responses
            end
        end
                
        % Picking mean number of trials from every participant in a group       
        mean_trnum = floor(mean(trnum_tot)); % computes the mean number of trials that needs to be picked up from each participant
       
        % Subject_Group.Combtr stores the concatenated trials (/ta/ and /pa/) of subjects in each group
        
  %{   
    Here we check for each participant in a group his/her combined (/ta/ and /pa/) trial number. 
    If the number is higher than the mean, we randomly pick the 'mean_trnum' number of trials from the subject. 
    This is concatenated and stored in Subject_Group.Combtr. 
    The trial index for the selected trials is stored in Subject_Group_Prestim.randnum
  %}
        
     for i = 1: Grpnum
         Subject_Group{i}.McG.Combtr = zeros(size(data{1}.McG_ta,1), size(data{1}.McG_ta, 2), 1);
         Subject_Group{i}.Cong_pa.tr = zeros(size(data{1}.McG_ta,1), size(data{1}.McG_ta, 2), 1);
         Subject_Group{i}.Cong_ta.tr = zeros(size(data{1}.McG_ta,1), size(data{1}.McG_ta, 2), 1);
         Subject_Group{i}.Cong_ka.tr = zeros(size(data{1}.McG_ta,1), size(data{1}.McG_ta, 2), 1);            
     end
        
        
     for grp = 1: Grpnum           
         Sub_indx = Subject_Group{grp}.Subindx;           
         for i = 1:length(Sub_indx)               
             if trnum_tot(Sub_indx(i)) <= mean_trnum 
                Subject_Group{grp}.McG.Combtr = cat(3,Subject_Group{grp}.McG.Combtr, data{Sub_indx(i)}.McG_ta, data{Sub_indx(i)}.McG_pa);
                Subject_Group{grp}.McG.randnum{i} = cell(1,1);
                Subject_Group{grp}.McG.Participant_tr{i} = cat(3, data{Sub_indx(i)}.McG_ta, data{Sub_indx(i)}.McG_pa);

                Subject_Group{grp}.Cong_pa.tr = cat(3, Subject_Group{grp}.Cong_pa.tr, data{Sub_indx(i)}.Cong_pa);
                Subject_Group{grp}.Cong_pa.randnum{i} = cell(1,1);
                Subject_Group{grp}.Cong_pa.Participant_tr{i} = data{Sub_indx(i)}.Cong_pa;

                Subject_Group{grp}.Cong_ta.tr = cat(3,Subject_Group{grp}.Cong_ta.tr, data{Sub_indx(i)}.Cong_ta);
                Subject_Group{grp}.Cong_ta.randnum{i} = cell(1,1);
                Subject_Group{grp}.Cong_ta.Participant_tr{i} = data{Sub_indx(i)}.Cong_ta;

                Subject_Group{grp}.Cong_ka.tr = cat(3,Subject_Group{grp}.Cong_ka.tr, data{Sub_indx(i)}.Cong_ka);
                Subject_Group{grp}.Cong_ka.randnum{i} = cell(1,1);
                Subject_Group{grp}.Cong_ka.Participant_tr{i} = data{Sub_indx(i)}.Cong_ka;
             end
                
             if trnum_tot(Sub_indx(i)) > mean_trnum
                 
                Comb_McG = cat(3, data{Sub_indx(i)}.McG_ta, data{Sub_indx(i)}.McG_pa);
                randnum_McG = randperm(size(Comb_McG,3));
                randnum_cong_pa = randperm(size(data{Sub_indx(i)}.Cong_pa,3));
                randnum_cong_ta = randperm(size(data{Sub_indx(i)}.Cong_ta,3));
                randnum_cong_ka = randperm(size(data{Sub_indx(i)}.Cong_ka,3));

                for tr = 1:mean_trnum
                    temp_McG(:,:,tr) = Comb_McG(:,:,randnum_McG(tr));
                end 

                % temp_Cong_pa
                if size(data{Sub_indx(i)}.Cong_pa,3) <= mean_trnum
                    temp_cong_pa = data{Sub_indx(i)}.Cong_pa;
                    Subject_Group{grp}.Cong_pa.randnum{i} = randnum_cong_pa;
                else
                    for tr = 1:mean_trnum
                        temp_cong_pa(:,:,tr) = data{Sub_indx(i)}.Cong_pa(:,:,randnum_cong_pa(tr));
                    end
                    Subject_Group{grp}.Cong_pa.randnum{i} = randnum_cong_pa(1:mean_trnum);
                end

            % temp_Cong_ta
                if size(data{Sub_indx(i)}.Cong_ta,3) <= mean_trnum
                    temp_cong_ta = data{Sub_indx(i)}.Cong_ta;
                    Subject_Group{grp}.Cong_ta.randnum{i} = randnum_cong_ta;
                else
                    for tr = 1:mean_trnum
                        temp_cong_ta(:,:,tr) = data{Sub_indx(i)}.Cong_ta(:,:,randnum_cong_ta(tr));
                    end
                    Subject_Group{grp}.Cong_ta.randnum{i} = randnum_cong_ta(1:mean_trnum);
                end

            % temp_Cong_ka
                if size(data{Sub_indx(i)}.Cong_ka,3) <= mean_trnum
                    temp_cong_ka = data{Sub_indx(i)}.Cong_ka;
                    Subject_Group{grp}.Cong_ka.randnum{i} = randnum_cong_ka;
                else
                    for tr = 1:mean_trnum
                        temp_cong_ka(:,:,tr) = data{Sub_indx(i)}.Cong_ka(:,:,randnum_cong_ka(tr));
                    end
                        Subject_Group{grp}.Cong_ka.randnum{i} = randnum_cong_ka;
                end
                Subject_Group{grp}.McG.Combtr  = cat(3,Subject_Group{grp}.McG.Combtr,temp_McG);
                Subject_Group{grp}.McG.randnum{i} = randnum_McG(1:mean_trnum);
                Subject_Group{grp}.McG.Participant_tr{i} = temp_McG;

                Subject_Group{grp}.Cong_pa.tr = cat(3,Subject_Group{grp}.Cong_pa.tr,temp_cong_pa);
                Subject_Group{grp}.Cong_pa.Participant_tr{i} = temp_cong_pa;

                Subject_Group{grp}.Cong_ta.tr = cat(3,Subject_Group{grp}.Cong_ta.tr,temp_cong_ta);
                Subject_Group{grp}.Cong_ta.Participant_tr{i} = temp_cong_ta;

                Subject_Group{grp}.Cong_ka.tr = cat(3,Subject_Group{grp}.Cong_ka.tr,temp_cong_ka);
                Subject_Group{grp}.Cong_ka.Participant_tr{i} = temp_cong_ka;                      
             end         
         end
     
            Subject_Group{grp}.McG.Combtr(:,:,1)=[];
            Subject_Group{grp}.Cong_pa.tr(:,:,1)=[];
            Subject_Group{grp}.Cong_ta.tr(:,:,1)=[];
            Subject_Group{grp}.Cong_ka.tr(:,:,1)=[];
            
            Trnum_Grp(grp) = size(Subject_Group{grp}.McG.Combtr,3);
            clear temp Comb randnum
      end
        
        
        % equalizing trials across groups
        fprintf('\nEqualizing trials across groups...')
        
        min_trnum_grp = min(Trnum_Grp);
        
        for i = 1:length(Subject_Group)
            if size(Subject_Group{i}.McG.Combtr,3) <= min_trnum_grp
                McG_grptrls_frGlobCoh{i}.tr = Subject_Group{i}.McG.Combtr;
                McG_grptrls_frGlobCoh{i}.rndnum = [];
                
                Cong_pa_grptrls_frGlobCoh{i}.tr = Subject_Group{i}.Cong_pa.tr;
                Cong_pa_grptrls_frGlobCoh{i}.rndnum = [];
                
                Cong_ta_grptrls_frGlobCoh{i}.tr = Subject_Group{i}.Cong_ta.tr;
                Cong_ta_grptrls_frGlobCoh{i}.rndnum = [];
                
                Cong_ka_grptrls_frGlobCoh{i}.tr = Subject_Group{i}.Cong_ka.tr;
                Cong_ka_grptrls_frGlobCoh{i}.rndnum = [];
            end
            if size(Subject_Group{i}.McG.Combtr,3) > min_trnum_grp
                McG_rnd_grp = randperm(size(Subject_Group{i}.McG.Combtr,3));
                cong_pa_rnd_grp = randperm(size(Subject_Group{i}.Cong_pa.tr,3));
                cong_ta_rnd_grp = randperm(size(Subject_Group{i}.Cong_ta.tr,3));
                cong_ka_rnd_grp = randperm(size(Subject_Group{i}.Cong_ka.tr,3));
                
                for tr = 1:min_trnum_grp
                    McG_grptrls_frGlobCoh{i}.tr(:,:,tr) = Subject_Group{i}.McG.Combtr(:,:,McG_rnd_grp(tr));
                    McG_grptrls_frGlobCoh{i}.rndnum = McG_rnd_grp(1:min_trnum_grp);
                    
                    Cong_pa_grptrls_frGlobCoh{i}.tr(:,:,tr) = Subject_Group{i}.Cong_pa.tr(:,:,cong_pa_rnd_grp(tr));
                    Cong_pa_grptrls_frGlobCoh{i}.rndnum = cong_pa_rnd_grp(1:min_trnum_grp);
                    
                    Cong_ta_grptrls_frGlobCoh{i}.tr(:,:,tr) = Subject_Group{i}.Cong_ta.tr(:,:,cong_ta_rnd_grp(tr));
                    Cong_ta_grptrls_frGlobCoh{i}.rndnum = cong_ta_rnd_grp(1:min_trnum_grp);
                    
                    Cong_ka_grptrls_frGlobCoh{i}.tr(:,:,tr) = Subject_Group{i}.Cong_ka.tr(:,:,cong_ka_rnd_grp(tr));
                    Cong_ka_grptrls_frGlobCoh{i}.rndnum = cong_ka_rnd_grp(1:min_trnum_grp);
                end
            end
        end
        
        fprintf('Done.\n');
        
        
        %{
        %%%% saving the group trials for cluster-based permutation testing for Global Coherence
        grptrls_frMnteCrlo.McG = McG_grptrls_frGlobCoh;
        grptrls_frMnteCrlo.Cong_pa = Cong_pa_grptrls_frGlobCoh;
        grptrls_frMnteCrlo.Cong_ta = Cong_ta_grptrls_frGlobCoh;
        grptrls_frMnteCrlo.Cong_ka = Cong_ka_grptrls_frGlobCoh;       
        
        save('C:\Users\dipanjanroy\Desktop\Project\Included\grptrls_InterSub_Coh_frMnteCrlo_pre800', 'grptrls_frMnteCrlo')       
        %}
        %McG_grptrls_frGlobCoh = grptrls_frMnteCrlo.McG;
        
        %%% Computing the Global Coherence combined groupwise
        % McGurk trials
        Grpnum = 2;
        Sc = cell(Grpnum,1); Cmat = cell(Grpnum,1); Ctot = cell(Grpnum,1);
        
        for i = 1:Grpnum
            fprintf('Computing Global Coherence of Prestimulus McGurk trials of group %d of %d...',i,Grpnum);
            
            [S,Cmt,Ct,~,Ce,f] = CrossSpecMatc(McG_grptrls_frGlobCoh{i}.tr, win, params);
            Sc{i}   = S;
            Cmat{i} = Cmt;
            Ctot{i} = Ct;
            Cent{i} = Ce;
            clear S Cmt Ct Ce
            fprintf('Done.\n');
        end
        
        GlobCoh.Group.McG.Sc = Sc;
        GlobCoh.Group.McG.Cmat = Cmat;
        GlobCoh.Group.McG.Ctot = Ctot;
        GlobCoh.Group.McG.Cent = Cent;
        GlobCoh.f = f;
        
        
        % Congruent /pa/ trials      
        Sc = cell(Grpnum,1); Cmat = cell(Grpnum,1); Ctot = cell(Grpnum,1);
        
        for i = 1:Grpnum
            fprintf('Computing Global Coherence of Prestimulus Congruent /pa/ trials of group %d of %d...',i,Grpnum);
            
            [S,Cmt,Ct,~,Ce,f] = CrossSpecMatc(Cong_pa_grptrls_frGlobCoh{i}.tr, win, params);
            Sc{i}   = S;
            Cmat{i} = Cmt;
            Ctot{i} = Ct;
            Cent{i} = Ce;
            clear S Cmt Ct Ce
            fprintf('Done.\n');
        end
        
        GlobCoh.Group.Cong_pa.Sc = Sc;
        GlobCoh.Group.Cong_pa.Cmat = Cmat;
        GlobCoh.Group.Cong_pa.Ctot = Ctot;
        GlobCoh.Group.Cong_pa.Cent = Cent;        
        GlobCoh.f = f;
        
        clear Ctot Cmat Sc Cent
        
        
        % Congruent /ta/ trials       
        Sc = cell(Grpnum,1); Cmat = cell(Grpnum,1); Ctot = cell(Grpnum,1);
        
        for i = 1:Grpnum
            fprintf('Computing Global Coherence of Prestimulus Congruent /ta/ trials of group %d of %d...',i,Grpnum);
            
            [S,Cmt,Ct,~,Ce,f] = CrossSpecMatc(Cong_ta_grptrls_frGlobCoh{i}.tr, win, params);
            Sc{i}   = S;
            Cmat{i} = Cmt;
            Ctot{i} = Ct;
            Cent{i} = Ce;
            clear S Cmt Ct Ce
            fprintf('Done.\n');
        end
        
        GlobCoh.Group.Cong_ta.Sc = Sc;
        GlobCoh.Group.Cong_ta.Cmat = Cmat;
        GlobCoh.Group.Cong_ta.Ctot = Ctot;
        GlobCoh.Group.Cong_ta.Cent = Cent;
        GlobCoh.f = f;
        
        clear Ctot Cmat Sc Cent
        
                
        % Congruent /ka/ trials       
        Sc = cell(Grpnum,1); Cmat = cell(Grpnum,1); Ctot = cell(Grpnum,1);
        
        for i = 1:Grpnum
            fprintf('Computing Global Coherence of Prestimulus Congruent /ka/ trials of group %d of %d...',i,Grpnum);
            
            [S,Cmt,Ct,~,Ce,f] = CrossSpecMatc(Cong_ka_grptrls_frGlobCoh{i}.tr, win, params);
            Sc{i}   = S;
            Cmat{i} = Cmt;
            Ctot{i} = Ct;
            Cent{i} = Ce;
            clear S Cmt Ct Ce
            fprintf('Done.\n');
        end
        
        GlobCoh.Group.Cong_ka.Sc = Sc;
        GlobCoh.Group.Cong_ka.Cmat = Cmat;
        GlobCoh.Group.Cong_ka.Ctot = Ctot;
        GlobCoh.Group.Cong_ka.Cent = Cent;
        GlobCoh.f = f;
        
        clear Ctot Cmat Sc Cent
        
        for grp = 1:Grpnum
            GlobCoh.trialinfo{grp} = Subject_Group{grp};            
        end
        
                
        %%% Participant-wise global coherence %%%
        
        [Sort_Percent_McG_ta,idx]= sort(Percent_McG_ta);
        for sub = 1: length(data)
            Percent_McG_ta(sub) = data{sub}.percent_McG.ta;
            Percent_McG_pa(sub) = data{sub}.percent_McG.pa;
        end
        
        X_ta_percnt = sort(Percent_McG_ta);
        X_ta_percnt = X_ta_percnt';
                
        Idx_from_grp = cat(1,Subject_Group{1}.Subindx, Subject_Group{2}.Subindx);
        
        Participantwise.McG.tr = cat(2, Subject_Group{1}.McG.Participant_tr, Subject_Group{2}.McG.Participant_tr);
        Participantwise.Cong_pa.tr = cat(2, Subject_Group{1}.Cong_pa.Participant_tr, Subject_Group{2}.Cong_pa.Participant_tr);
        Participantwise.Cong_ta.tr = cat(2, Subject_Group{1}.Cong_ta.Participant_tr, Subject_Group{2}.Cong_ta.Participant_tr);
        Participantwise.Cong_ka.tr = cat(2, Subject_Group{1}.Cong_ka.Participant_tr, Subject_Group{2}.Cong_ka.Participant_tr);
        
        
        for i = 1:length(data)
            x = find(Idx_from_grp == idx(i));
            McG_tr{i} = Participantwise.McG.tr{x};
            Cong_pa_tr{i} = Participantwise.Cong_pa.tr{x};
            Cong_ta_tr{i} = Participantwise.Cong_ta.tr{x};
            Cong_ka_tr{i} = Participantwise.Cong_ka.tr{x};
        end
                
        clear Participantwise
        
        % Subjectwise McG
        for i = 1:length(data)
            fprintf('Computing Global-Coherence of Prestimulus McGurk trials of participant %d of %d...',i,length(data));
            [~,~,Ctot,~,~,f] = CrossSpecMatc(McG_tr{i},win,params);
            Data_sort{i}.McG.Ctot = Ctot;
            Data_sort{i}.McG.tr   = McG_tr{i};
            fprintf('Done.\n');
        end
        clear McG_tr
        
        % Subjectwise Cong/pa/
        for i = 1:length(data)
            fprintf('Computing Global-Coherence of Prestimulus Congruent /pa/ trials of participant %d of %d...',i,length(data));
            [~,~,Ctot,~,~,f] = CrossSpecMatc(Cong_pa_tr{i},win,params);
            Data_sort{i}.Cong_pa.Ctot = Ctot;
            Data_sort{i}.Cong_pa.tr   = Cong_pa_tr{i};
            fprintf('Done.\n');
        end
        clear Cong_pa_tr
        
        % Subjectwise Cong/ta/
        for i = 1:length(data)
            fprintf('Computing Global-Coherence of Prestimulus Congruent /ta/ trials of participant %d of %d...',i,length(data));
            [~,~,Ctot,~,~,f] = CrossSpecMatc(Cong_ta_tr{i},win,params);
            Data_sort{i}.Cong_ta.Ctot = Ctot;
            Data_sort{i}.Cong_ta.tr   = Cong_ta_tr{i};
            fprintf('Done.\n');
        end
        clear Cong_ta_tr
        
        % Subjectwise Cong/ka/
        for i = 1:length(data)
            fprintf('Computing Global-Coherence of Prestimulus Congruent /ka/ trials of participant %d of %d...',i,length(data));
            [~,~,Ctot,~,~,f] = CrossSpecMatc(Cong_ka_tr{i},win,params);
            Data_sort{i}.Cong_ka.Ctot = Ctot;
            Data_sort{i}.Cong_ka.tr   = Cong_ka_tr{i};
            fprintf('Done.\n');
        end
        
        clear Cong_ka_tr
        
        GlobCoh.Participantwise = Data_sort;
                        
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
           
    case 'InterTrial_Groupwise'

        % The choice computes the global coherence of the percepts /ta/ and /pa/ within 
        % groups (grouped based on Participantwise McGurk percept)
        
        for i = 1:length(data)
            Percent_McG_ta(i)= data{i}.percent_McG.ta;
        end
        
        [Sort_Percent_McG_ta,idx] = sort(Percent_McG_ta);
        
        grpnum = 100/cfg.interval;
        
        grps = cell(grpnum,1);
        
        % computes the number of participants within a group based on the
        % specified susceptibility threshold in cfg.interval
        strt_percent = 0;
        end_percent = cfg.interval;
        
        for i = 1:grpnum
            if end_percent < 100
                grps{i} = find(Percent_McG_ta >= strt_percent & Percent_McG_ta < end_percent);
            else
                grps{i} = find(Percent_McG_ta >= strt_percent & Percent_McG_ta <= end_percent);
            end
            strt_percent = end_percent;
            end_percent = end_percent + cfg.interval;
        end
        
        %Extracting the Behavior(/ta/ and /pa/ %ages) and trials from each subject
        sort_data = data(idx);
        
        Sub_len = length(sort_data);       % number of subjects
        Percent_McG_ta = zeros(Sub_len,1);    % stores the %age of /ta/ response of each subject
        Percent_McG_pa = zeros(Sub_len,1);    % stores the %age of /pa/ response of each subject
        trnum_ta  = zeros(Sub_len,1);         % stores the NUMBER of /ta/ trials of each subject
        trnum_pa  = zeros(Sub_len,1);         % stores the NUMBER of /pa/ trials of each subject
        trnum_tot = zeros(Sub_len,1);         % stores the total NUMBER of McG trials (/ta/+/pa/)  of each subject
        
        for sub = 1 : length(sort_data)
            Percent_McG_ta(sub) = sort_data{sub}.percent_McG.ta;
            Percent_McG_pa(sub) = sort_data{sub}.percent_McG.pa;
            trnum_ta(sub) = size(sort_data{sub}.McG_ta,3);
            trnum_pa(sub) = size(sort_data{sub}.McG_pa,3);
        end
        
        trnum_tot = trnum_ta + trnum_pa; % array that stores the total number of trials(ta+pa) of each subject
        
        
        % variables to store the concatenated /pa/ and /ta/ trials based on groups
        
        for i = 1:grpnum
           ta_tr{i} = zeros(801,64,1);
           pa_tr{i} = zeros(801,64,1);
           cong_pa_tr{i} = zeros(801,64,1);
           cong_ta_tr{i} = zeros(801,64,1);
           cong_ka_tr{i} = zeros(801,64,1);    
        end
        
        % Concatenating the trials
        
        for i = 1:grpnum
            for j = 1:length(grps{i})
                ta_tr{i} = cat(3,ta_tr{i}, sort_data{grps{i}(j)}.McG_ta);
                pa_tr{i} = cat(3,pa_tr{i}, sort_data{grps{i}(j)}.McG_pa);
                cong_pa_tr{i} = cat(3,cong_pa_tr{i}, sort_data{grps{i}(j)}.Cong_pa);
                cong_ta_tr{i} = cat(3,cong_ta_tr{i}, sort_data{grps{i}(j)}.Cong_ta);
                cong_ka_tr{i} = cat(3,cong_ka_tr{i}, sort_data{grps{i}(j)}.Cong_ka);
            end
            ta_tr{i}(:,:,1) = [];
            pa_tr{i}(:,:,1) = []; 
            cong_pa_tr{i}(:,:,1) = [];
            cong_ta_tr{i}(:,:,1) = [];
            cong_ka_tr{i}(:,:,1) = [];
        end
        
        % finding the minimum number of trials for randomly picking from each condition
        % the McG and congruent trial numbers are ordered as /pa/,/ta/,/cong_pa/,/cong_ta/,/cong_ka/...
        z = 1;
        for i = 1:grpnum
            trnum(z) = size(pa_tr{i},3);
            trnum(z+1)= size(ta_tr{i},3);
            trnum(z+2)= size(cong_pa_tr{i},3);
            trnum(z+3)= size(cong_ta_tr{i},3);
            trnum(z+4)= size(cong_ka_tr{i},3);
            z=z+5;
        end
        
        % minimum number of trials to be selected from each condition
        min_trnum = min(trnum);
        
        % McGurk /pa/,/ta/ and Cong /pa/,/ta/,/ka/ are order in trnum as [/pa/,/ta/,/cong_pa/,/cong_ta/,/cong_ka/...]
        z = 1;
        for i = 1:grpnum
            pa_rnd{i} = randperm(trnum(z));
            ta_rnd{i} = randperm(trnum(z+1));
            cong_pa_rnd{i} = randperm(trnum(z+2));
            cong_ta_rnd{i} = randperm(trnum(z+3));
            cong_ka_rnd{i} = randperm(trnum(z+4));
            z = z+5;
        end
     
        % Randomly picking trials from the concatenated trials from each group
        
        fprintf('Randomly picking %d trials from each condition...', min_trnum);
        
        for i = 1 : grpnum
            for j = 1:min_trnum
                pa_tr_rnd{i}(:,:,j) = pa_tr{i}(:,:,pa_rnd{i}(j));
                ta_tr_rnd{i}(:,:,j) = ta_tr{i}(:,:,ta_rnd{i}(j));
                cong_pa_tr_rnd{i}(:,:,j) = cong_pa_tr{i}(:,:,cong_pa_rnd{i}(j));
                cong_ta_tr_rnd{i}(:,:,j) = cong_ta_tr{i}(:,:,cong_ta_rnd{i}(j));
                cong_ka_tr_rnd{i}(:,:,j) = cong_ka_tr{i}(:,:,cong_ka_rnd{i}(j));
            end
        end
        
        fprintf('Done.\n');
        
        %{
          % saving the group trials for cluster-based permutation testing for Global Coherence
            grptrls_frMnteCrlo.McG.pa = pa_tr_rnd;
            grptrls_frMnteCrlo.McG.ta = ta_tr_rnd;
            grptrls_frMnteCrlo.Cong_pa = cong_pa_tr_rnd;
            grptrls_frMnteCrlo.Cong_ta = cong_ta_tr_rnd;
            grptrls_frMnteCrlo.Cong_ka = cong_ka_tr_rnd;
        
           save('C:\Users\dipanjanroy\Desktop\Project\Included\grptrls_InterTrial_Coh_frMnteCrlo_pre800', 'grptrls_frMnteCrlo') 
           %save('C:\Users\dipanjanroy\Desktop\Project\Included\grptrls_InterTrial_Coh_frMnteCrlo_pre500', 'grptrls_frMnteCrlo')
              
        %}
        
        %preparing the output structure
        
        for i = 1:grpnum
            Group{i}.McG.pa.trials = pa_tr_rnd{i};
            Group{i}.McG.pa.trnum  = pa_rnd{i};
            
            Group{i}.McG.ta.trials = ta_tr_rnd{i};
            Group{i}.McG.ta.trnum  = ta_rnd{i};
            
            Group{i}.Cong_pa.trials = cong_pa_tr_rnd{i};
            Group{i}.Cong_pa.trnum = cong_pa_rnd{i};
            
            Group{i}.Cong_ta.trials = cong_ta_tr_rnd{i};
            Group{i}.Cong_ta.trnum = cong_ta_rnd{i};
            
            Group{i}.Cong_ka.trials = cong_ka_tr_rnd{i};
            Group{i}.Cong_ka.trnum = cong_ka_rnd{i};           
        end
        
        
        % Computing the global coherence        
        for i = 1:grpnum
            fprintf('Computing Global Coherence of Prestimulus McG /pa/ trials of group %d of %d...',i,grpnum);
            [S,Cmt,Ct,~,~,f] = CrossSpecMatc(Group{i}.McG.pa.trials,win,params);
            Group{i}.McG.pa.Sc   = S;
            Group{i}.McG.pa.Cmat = Cmt;
            Group{i}.McG.pa.Ctot = Ct;           
            clear S Cmt Ct
            fprintf('Done.\n');
        end
        
        
        for i = 1:grpnum
            fprintf('Computing Global Coherence of Prestimulus McG /ta/ trials of group %d of %d...',i,grpnum);
            [S,Cmt,Ct,~,~,f] = CrossSpecMatc(Group{i}.McG.ta.trials,win,params);
            Group{i}.McG.ta.Sc   = S;
            Group{i}.McG.ta.Cmat = Cmt;
            Group{i}.McG.ta.Ctot = Ct;
            Group{i}.f = f;
            clear S Cmt Ct
            fprintf('Done.\n');
        end
        
        
        for i = 1:grpnum
            fprintf('Computing Global Coherence of Prestimulus Congruent /pa/ trials of group %d of %d...',i,grpnum);
            [S,Cmt,Ct,~,~,f] = CrossSpecMatc(Group{i}.Cong_pa.trials,win,params);
            Group{i}.Cong_pa.Sc   = S;
            Group{i}.Cong_pa.Cmat = Cmt;
            Group{i}.Cong_pa.Ctot = Ct;
            clear S Cmt Ct
            fprintf('Done.\n');
        end
        
        for i = 1:grpnum
            fprintf('Computing Global Coherence of Prestimulus Congruent /ta/ trials of group %d of %d...',i,grpnum);
            [S,Cmt,Ct,~,~,f] = CrossSpecMatc(Group{i}.Cong_ta.trials,win,params);
            Group{i}.Cong_ta.Sc   = S;
            Group{i}.Cong_ta.Cmat = Cmt;
            Group{i}.Cong_ta.Ctot = Ct;
            clear S Cmt Ct
            fprintf('Done.\n');
        end
        
        for i = 1:grpnum
            fprintf('Computing Global Coherence of Prestimulus Congruent /ka/ trials of group %d of %d...',i,grpnum);
            [S,Cmt,Ct,~,~,f] = CrossSpecMatc(Group{i}.Cong_ka.trials,win,params);
            Group{i}.Cong_ka.Sc   = S;
            Group{i}.Cong_ka.Cmat = Cmt;
            Group{i}.Cong_ka.Ctot = Ct;
            clear S Cmt Ct
            fprintf('Done.\n');
        end
        
        GlobCoh = Group;

end

end  % end of the function

