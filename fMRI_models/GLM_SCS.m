function GLM_SCS(sub, cononly)
%fMRI dietary self-control model to compare to Hare 2009 original analysis

% Inputs: 
% sub = participant identifier as a character string (e.g. '9999')
% cononly = enter 1: run only the contrasts, 0 = estimate the SPM

% abbreviations: SC = self-control

% acquisition time (1 volume)
TR = 2.344;

% specifics for this task:
nrTrials = 100; % total number of trials in the paradigm
ndyns = 354; % number of volumes acquired
% task prefix ('fc'= food choice, 'er' = emotion regulation)
task_prefix = 'fc';

% name of the model
model_name = 'GLM_SCS';

% orthogonalize or not?
orth_switch = 0; % 0 = do not orthogonalize pmods, 1 = orthogonalize pmods

% include onset only if it has at least 3 events, otherwise estimation becomes too noisy
min_event_num = 3; % want to have at least one entry for onset + number of pmods on each regressor

% paths
image_folder = '/Users/username/mrdata/ESC/subs/SNS_ESC1_'; % folder with preprocessed functional images
behavior_path = '/Users/username/mrdata/ESC/behavior/'; % all behavioral files of all participants live in this folder
physio_path = '/Users/username/mrdata/ESC/hrdata/'; % this is where the physio files live
save_path = ['/Users/username/mrdata/ESC/results/' model_name]; % where model results are saved
save_path2 = [save_path filesep sub]; % folder to save each participant's results
% complete path to folder in which preprocessed images live for this task
H_folder = [image_folder sub filesep 'functional' filesep task_prefix filesep];

% create and fill onsets and parametric modulators
nsessions = 1; % number of runs


if cononly == 0 % want to estimate the SPM (== 1 runs contrasts only)
    
    %% define design conditions
    
    mov_reg{nsessions} = {};
    filestr{nsessions} = {};
    ons_duration{nsessions} = {};
    ons{nsessions} = {};
    ons_modulate{nsessions} = {};
    
    % onsets
    
    ons_name{1} = 'SCsuccess'; % participant used self-control
    ons_name{2} = 'SCfailure'; % participant did not use self-control
    ons_name{3} = 'NoChall'; % this includes no challenge trials and trials in which an item was neutral on either taste or health or both
    ons_name{4} = 'Missed'; % no response given on this trial
   
    
    % get values
    imglist = load([behavior_path 'imagelist_' sub '.mat']);
    alltaste = load([behavior_path 'results_taste_' sub '.mat'],'responses');
    allhealth = load([behavior_path 'results_health_' sub '.mat'],'responses');
       
    for session = 1:nsessions
        
        
        % create vectors of behavioral data
        ss = load([behavior_path 'food_choices_' sub '.mat']);
        nonmiss = ~isnan(ss.responses); % pictures for which an answer was given
        miss = isnan(ss.responses); % missed trial
        
        
        tr = zeros(nrTrials,1); % taste ratings for the foods
        hr = zeros(nrTrials,1); % health ratings for the foods
        scs = zeros(nrTrials,1); % identifier for SC success trials
        scf = zeros(nrTrials,1); % identifier for SC failure trials
        noc = zeros(nrTrials,1); % identifier for no challenge trials
        rts = zeros(nrTrials,1); % reaction time on that trial
        
        %  read in ratings - ratings only work for non-missed trials
        for i = 1 : nrTrials
            
            % construct vectors with values for modulators
            if ~isnan(ss.responses(i))
                tr(i) = alltaste.responses(imglist.imgvect.blocks(i));
                hr(i) = allhealth.responses(imglist.imgvect.blocks(i));
                       
                % neutral lower and upper bounds around 0 as they were
                % depicted on the rating scale (neutral area are 5 percent on each side
                % of zero); rating scale raned from 0 to 426 pixels 
                neutralbottom = 213-21;
                neutraltop = 213+21;
                
                % simple yes/no classification
                if ss.responses(i) == 1
                    yes(i) = 1;
                    no(i) = 0;
                end
                if ss.responses(i) == 2
                    yes(i) = 0;
                    no(i) = 1;
                end
                
                % classify by self-control category
                
                % Challenges
                
                % tasty-unhealthy case
                % item is tastier and less healthy and chosen
                if (tr(i) > neutraltop) && (hr(i) < neutralbottom) && (ss.responses(i) == 1)
                    scf(i) = 1; % failure
                end
                
                % item is tastier and less healthy and not chosen
                if (tr(i) > neutraltop) && (hr(i) < neutralbottom) && (ss.responses(i) == 2)
                    scs(i) = 1; %success
                end
                
                % healthy-untasty case
                % item is less tasty and healthier and chosen
                if (tr(i) < neutralbottom) && (hr(i) > neutraltop) && (ss.responses(i) == 1)
                    scs(i) = 1; %success
                end
                
                %item is less tasty and healthier and not chosen
                if (tr(i) < neutralbottom) && (hr(i) > neutraltop) && (ss.responses(i) == 2)
                    scf(i) = 1; %failure
                end
                
                % No Challenges
                
                % healthy-tasty case
                % item is tastier and healthier and chosen
                if (tr(i) > neutraltop) && (hr(i) > neutraltop) && (ss.responses(i) == 1)
                    noc(i) = 1; %no challenge
                end
                
                % item is tastier and healthier and not chosen
                if (tr(i) > neutraltop) && (hr(i) > neutraltop) && (ss.responses(i) == 2)
                    noc(i) = 1; %no challenge
                end
                
                % unhealthy-untasty case
                % item is less tasty and less healthy and chosen
                if (tr(i) < neutralbottom) && (hr(i) < neutralbottom) && (ss.responses(i) == 1)
                    noc(i) = 1; %no challenge
                end
                
                % item is less tasty and less healthy and not chosen
                if (tr(i) < neutralbottom) && (hr(i) < neutralbottom) && (ss.responses(i) == 2)
                    noc(i) = 1; %no challenge
                end
                
                % Neutral Item goes into No Challenge count
                
                % taste in the neutral zone
                if (tr(i) > neutralbottom) && (tr(i) < neutraltop)
                    noc(i) = 1; %neutral -> goes into no challange
                end
                
                % health in the neutral zone
                if (hr(i) > neutralbottom) && (hr(i) < neutraltop)
                    noc(i) = 1; %neutral -> goes into no challenge
                end
                
            end
            
        end %end trial loop to populate pmod values
        
        
        %% set up regressors of interest
        
        % subtract the absolute start_time (trigger)
        % from all time stamps, so the time line of the onsets starts at 0
        ss.timing.trial_start_times = ss.timing.trial_start_times-ss.timing.start_time; % start time for each choice
        
        % onsets
        ons{1,session} = ss.timing.trial_start_times(scs > 0 & nonmiss); % all SC success
        ons{2,session} = ss.timing.trial_start_times(scf > 0 & nonmiss); % all SC failure
        ons{3,session} = ss.timing.trial_start_times(noc > 0 & nonmiss); % all no challenge
        ons{4,session} = ss.timing.trial_start_times(miss); % all missed trials
        
        % durations  - as vectors
        ons_duration{1,session} = ss.timing.reaction_times(scs > 0 & nonmiss)'; % duration equals reaction time
        ons_duration{2,session} = ss.timing.reaction_times(scf > 0 & nonmiss)'; % duration equals reaction time
        ons_duration{3,session} = ss.timing.reaction_times(noc > 0 & nonmiss)'; % duration equals reaction time
        maxtrialtimes = 3*ones(nrTrials,1); % missed trials are modeled as the full 3 sec that participants had to decide
        ons_duration{4,session} = maxtrialtimes(miss)'; 
        

        
        %% movement - find bad onsets and define separate regressor for 3 TRs before and 1 after
        
        clear mov_reg{session}; % careful! do NOT clear mov_reg as a whole or will remove needed data for session > 1
        
        rpfile=dir([H_folder '/rp_*.txt']);
        mov_reg{session} = load([H_folder rpfile(session).name]);
        
        % model out movement that exceeds 2 mm in the rotations or 2 dg in the translations
        % create regressor to account for anything that is related
        % to this movement (over and above the onsets/events you have modeled)
        
        % careful! this is assumes that your movement modeling has not been done with the physIO toolbox!
        pos_excess_mov = mov_reg{session} > 2;
        neg_excess_mov = mov_reg{session} < -2;
        total_excess_mov = pos_excess_mov | neg_excess_mov;
        
        % test if any of the rows of the excess movement regressor contains an entry
        stick_regressor = any(total_excess_mov,2);
        
        % index and label excessive movement events (incl 3 before and 1
        % after)
        xcsmovindex = find(stick_regressor); % get the position/onsets of all excess movement events
        
        if length(xcsmovindex) >= 1 % just do this if there is any excess movement
            for i = 1:length(xcsmovindex)
                if xcsmovindex(i) == 1 % if this happens on the first trial
                    stick_regressor(xcsmovindex(i):xcsmovindex(i)+1)=1;
                elseif xcsmovindex(i) == 2 % ... second trial
                    stick_regressor(xcsmovindex(i)-1:xcsmovindex(i)+1)=1;
                elseif xcsmovindex(i) == 3 % ... or third trial ...
                    stick_regressor(xcsmovindex(i)-2:xcsmovindex(i)+1)=1;
                else % always go 3 trials back and flag them
                    stick_regressor(xcsmovindex(i)-3:xcsmovindex(i)+1)=1;
                end
            end
        end
        
        
        %% make physio and motion nuisance regressor (multiple regressor) files
        
        try
            % output of retroicor
            tmp = load([physio_path filesep 'SNS_ESC1_'  sub filesep task_prefix filesep 'physio_output' filesep task_prefix '_multiple_regressors.mat']);
            physio{session} = tmp.R;
            
            % put together physio and movement regressors
            if any(stick_regressor)
                mov_reg{session} = cat(2, physio{session}((1:ndyns),:), mov_reg{session}(1:ndyns,:), stick_regressor(1:ndyns) );
            else
                mov_reg{session} = cat(2, physio{session}((1:ndyns),:), mov_reg{session}(1:ndyns,:));
            end
            
        catch
            disp('Problem with retroicor for this participant, using only motion parms')
            
            retroicor = 0;
            if any(stick_regressor)
                mov_reg{session} = cat(2, mov_reg{session}(1:ndyns,:), stick_regressor(1:ndyns));
            else
                mov_reg{session} = mov_reg{session}(1:ndyns,:);
            end
            
        end
        
        R = mov_reg{session};
        
    end
    save([physio_path filesep 'SNS_ESC1_'  sub filesep task_prefix filesep 'physio_output' filesep task_prefix '_multiple_regressors_new.mat'], 'R')
    
    % complete path to multiple_regressors.mat file
    reg_file_name = fullfile([physio_path filesep 'SNS_ESC1_'  sub filesep task_prefix filesep 'physio_output' filesep ], [task_prefix '_multiple_regressors_new.mat']);
    
    %% make multiple conditions file
    
    names = ons_name; % names of conditions
    onsets = ons; % as defined above, all onsets
    durations = ons_duration; % duration = reaction time for each event
    
    % fill information for the parametric modulators for each onset
    pmod(1).name = {''}; % no pmod
    pmod(1).param = {[]}; % actual values of the parametric modulator
    pmod(1).poly = {1}; % polynomial expansions, keep on 1 unless you want also quadratic effects
    
    pmod(2).name = {''}; % no pmod
    pmod(2).param = {[]};
    pmod(2).poly = {1};
    
    pmod(3).name = {''};% no pmod
    pmod(3).param = {[]};
    pmod(3).poly = {1};
    
    pmod(4).name = {''}; % no pmod, nuisance regressor
    pmod(4).param = {[]};
    pmod(4).poly = {1};
    
    
    
    % save multiple_condition_file.mat
    fname2 = [save_path filesep 'MCF'];
    mkdir(fname2);
    
    fname = [save_path filesep 'MCF' filesep sub '_MCF_' task_prefix '.mat'];
    save (fname, 'names', 'onsets' ,'durations' ,'pmod') ;

    
    %% Module 1 for the SPM batch: populate the SPM design
    
    m = 1;
    
    % check which onsets are too sparsely populated to estimate properly
    % inspect conditions: are there enough events?
    t = 1;
    t2 = 1;
    temp5 = 1;
    low_event_matrix = [];
    low_event_matrix_pmod = [];
    
    for p = 1:length(names)
        % only if the condition has more then the minimum number of events
        if length(onsets{p}) >= min_event_num 
            matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).name = names{p};
            matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).onset = onsets{p};
            matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).duration = durations{p};
            matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).tmod = 0;
            
            % orthogonalize pmods or not? 0 = no, 1 = yes
            matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).orth = orth_switch;
            
            % pmods
            z = 1;
            for w = 1:length(pmod(p).name)
                % only include pmods that had at least minimum number of events
                if ((length(unique(pmod(p).param{1, w})) >= min_event_num) == 1)
                    matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).pmod(z).name = pmod(p).name{w};
                    matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).pmod(z).param = pmod(p).param{w};
                    matlabbatch{m}.spm.stats.fmri_spec.sess(1).cond(t).pmod(z).poly = pmod(p).poly{w};
                    z = z+1;
                else
                    % write low event matrix - to know if you have too few events in this condition
                    low_event_matrix_pmod{t2,1} = sub;
                    low_event_matrix_pmod{t2,2} = pmod(p).name{w};
                    low_event_matrix_pmod{t2,3} = length(pmod(p).param{w});
                    low_event_matrix_pmod{t2,4} = pmod(p).param{w};
                    t2 = t2+1;
                end
            end
            t=t+1;
        else
            % write out
            low_event_matrix{temp5,1} = sub;
            low_event_matrix{temp5,2} = names{p};
            low_event_matrix{temp5,3} = length(onsets{p});
            temp5 = temp5+1;
        end
    end;
    
    
    % general parameters / SPM defaults for the first level model
    matlabbatch{m}.spm.stats.fmri_spec.dir = {save_path2};
    matlabbatch{m}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{m}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{m}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{m}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{m}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{m}.spm.stats.fmri_spec.sess(1).hpf = 128;
    matlabbatch{m}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{m}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{m}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{m}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{m}.spm.stats.fmri_spec.mthresh = 0.8;
    %run within explicit mask of skull-striped brain:
    matlabbatch{m}.spm.stats.fmri_spec.mask = {'/Applications/spm12/canonical/single_subj_T1_NoSkull_binary.nii'};
    matlabbatch{m}.spm.stats.fmri_spec.cvi = 'AR(1)'; % autocorrelation matrix
    
    % load scans
    % scans -> find all volumes and put them into structure
    scans = strcat(H_folder,cellstr(spm_select('ExtList',H_folder,'swausn_*',inf)));
    number_of_TRs = ndyns; %length(scans)
    matlabbatch{m}.spm.stats.fmri_spec.sess(1).scans = scans;
    
    % load physio
    matlabbatch{m}.spm.stats.fmri_spec.sess(1).multi_reg = {reg_file_name};

    
    %% Module 2 for the SPM batch: save and estimate the SPM
    
    m = m+1;
    
    spm_jobman('initcfg');
    save_path2=[save_path filesep sub];
    
    matlabbatch{m}.spm.stats.fmri_est.spmmat = {[save_path2 filesep 'SPM.mat']};
    matlabbatch{m}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{m}.spm.stats.fmri_est.method.Classical = 1; % use Restricted Maximum Likelihood method
    
    spm_jobman('run',matlabbatch)
    
    clear matlabbatch
    
end

%% First level contrasts (start over with a new batch, start again at module 1)

m = 1;

% contrasts
CON_NAMES = { 'SCsuccess>SCfailure',    'SCsuccess>NoChall',   ...
    };

% These two lists must have the same length. If you have 1 sample ttests, then use something like
% null for the negative conname. (can be any name that will not match one of the regressors)
connames_pos = { 'SCsuccess*bf',       'SCsuccess*bf',   ...
    };
connames_neg = { 'SCfailure*bf',     'NoChall*bf',   ...
    };


for j = 1:length(CON_NAMES)
    
    CONS = create_cons(connames_pos, connames_neg, [save_path2 filesep 'SPM.mat']);
end

sumcons = sum(CONS,2);
goodcons = ~isnan(sumcons);

matlabbatch{m}.spm.stats.con.spmmat = {[save_path2 filesep 'SPM.mat']};

CON_NAMES = CON_NAMES(goodcons);
CONS = CONS(goodcons,:);

for i = 1:length(CON_NAMES)
    matlabbatch{m}.spm.stats.con.consess{i}.tcon.name = CON_NAMES{i};
    matlabbatch{m}.spm.stats.con.consess{i}.tcon.convec = CONS(i,:);
    matlabbatch{m}.spm.stats.con.consess{i}.tcon.sessrep = 'none';  % don't replicate over sessions (replsc) if using create_cons function b/c already done
end
matlabbatch{m}.spm.stats.con.delete = 1; % this clears out all existing contrasts if the contrasts are re-run


% execute job and create con.nii files
spm_jobman('run', matlabbatch)

save([save_path2 filesep 'con_names_' sub], 'CON_NAMES') % save images for re-labeling contrasts later if cons get dropped bc not populated
save model_cons CON_NAMES; % save overview of all cons in the model results folder

    function con = create_cons(connames_pos, connames_neg, SPM)  % adapted from Jan Glascher's code (I think)
        
        temp = load(SPM);
        SPM = temp.SPM;
        
        conweight = 1;
        negweight = -1;
        
        con = zeros(length(connames_pos),length(SPM.xX.name));
        negcon = zeros(length(connames_pos),length(SPM.xX.name));
        
        for c = 1:length(connames_pos)
            
            for r = 1:length(SPM.xX.name)
                if strfind(SPM.xX.name{r}, connames_pos{c})
                    con(c,r) = conweight;
                end
            end
            
            
            for r = 1:length(SPM.xX.name)
                if strfind(SPM.xX.name{r}, connames_neg{c})
                    negcon(c,r) = negweight;
                end
            end
            
            
            % Normalize contrast weights to account for inequalities across sessions
            con(c,:) = con(c,:) ./ sum(con(c,:));
            negcon(c,:) = -1*(negcon(c,:) ./ min(-1,sum(negcon(c,:)))); %Multiplying by -1 preserves the sign, min(-1,sum) prevents division by 0
            con(c,:) = con(c,:) + negcon(c,:);
            
        end % end con weights loop
        
    end % end con function

end