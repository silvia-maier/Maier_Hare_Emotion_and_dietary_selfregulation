function GLM_ER(sub, cononly)
% fMRI first-level analysis for the main emotion regulation model

% Inputs: 
% sub = participant identifier as a character string (e.g. '9999')
% cononly = enter 1: run only the contrasts, 0 = estimate the SPM


% acquisition time (1 volume)
TR = 2.344;

% specifics for this task:
nrTrials = 100; % number of trials in the paradigm
ndyns = 681; % number of volumes acquired
if strcmp(sub, '1003') 
    ndyns = 679; % this participant was measured before we agreed to measure 681 volumes to have a margin
end

% task prefix ('fc'= food choice, 'er' = emotion regulation)
task_prefix = 'er';

% name of the model
model_name = 'GLM_ER';

% orthogonalize or not?
orth_switch = 0; % 0 = do not orthogonalize pmods, 1 = orthogonalize pmods

% include onset only if it has at least 3 events, otherwise estimation becomes too noisy
min_event_num = 3; % want to have at least one entry for onset + number of pmods on each regressor

% paths
image_folder = '/Users/username/mrdata/ESC/subs/SNS_ESC1_'; % location of the preprocessed functional images
behavior_path = '/Users/username/mrdata/ESC/behavior/'; % location of the behavioral data
physio_path = '/Users/username/mrdata/ESC/hrdata/'; % location of the physio files that contain ECG and breathing recordings
save_path = ['/Users/username/mrdata/ESC/results/' model_name]; % where model results are saved
save_path2 = [save_path filesep sub]; % folder to save each participant's results
% complete path to folder in which preprocessed images live for this task
H_folder = [image_folder sub filesep 'functional' filesep task_prefix filesep];

% create and fill onsets and parametric modulators
nsessions = 1; % number of runs

if cononly == 0 % if want to estimate the SPM, otherwise == 1 estimates only contrasts
    %% First-level model setup
    
    % define design conditions
    
    mov_reg{nsessions} = {};
    filestr{nsessions} = {};
    ons_duration{nsessions} = {};
    ons{nsessions} = {};
    ons_modulate{nsessions} = {};
    
    % onsets
    ons_name{1} = 'posViewEmo'; % positive picture viewing
    ons_name{2} = 'posRegSuccess'; % positive picture regulation success
    ons_name{3} = 'posRegFail'; % positive picture regulation failure
    ons_name{4} = 'negViewEmo'; % negative picture viewing
    ons_name{5} = 'negRegSuccess'; % negative picture regulation success
    ons_name{6} = 'negRegFail'; % negative picture regulation failure
    ons_name{7} = 'neuView'; % neutral picture viewing
    ons_name{8} = 'ratingScreen'; % participants rated their current mood
    
    % get values
    imglist = load([behavior_path 'imagelist_er_' sub '.mat']); % imgvect.blocks gives the picture number of the regulated images
    
    for session = 1:nsessions
        
        % create vectors of behavioral data
        ss = load([behavior_path 'er_scan_' sub '.mat']); % load rating data
        nonmiss = ~isnan(ss.responses.responses); % pictures for which a rating was made
        
        psr = load([behavior_path 'er_postscan_ratings_' sub '.mat']); % postscan ("view") ratings for the regulated images
        post_indexed_ratings = cat(2,psr.image_order, psr.responses.responses); % sort in order of the presentation to calculate regulation success
        
        post_rating_reg = NaN(nrTrials,1);
        % create a vector of sorted post scan ratings that match the
        % order presented in the scanner
        for i = 1:length(post_indexed_ratings)
            for j = 1:nrTrials
                if post_indexed_ratings(i,1) == imglist.imgvect.blocks(j,1)
                    post_rating_reg(j) = post_indexed_ratings(i,2);
                end
            end
        end
        
        er = zeros(nrTrials,1); % emotion ratings for the presented IAPS pictures
        rts = zeros(nrTrials,1); % reaction times for emotion ratings
        
        posViewEmo = zeros(nrTrials,1); % identifier for positive view trials
        posRegSuccess = zeros(nrTrials,1); % identifier for positive regulate success trials
        posRegFail = zeros(nrTrials,1); % identifier for positive regulate failure trials
        negViewEmo = zeros(nrTrials,1); % identifier for negative view trials
        negRegSuccess = zeros(nrTrials,1); % identifier for negative regulate success trials
        negRegFail = zeros(nrTrials,1); % identifier for negative regulate failure trials
        neuView = zeros(nrTrials,1); % identifier for neutral view trials
        
        % identify blocks of pos / neg pre-allocated trials:
        % we used 5 pseudo-random orders
        
        % for all success indices: initialize as one and only set to zero if fail
        % for all fail indices: initialize as one and set to zero if success
        
        switch imglist.rand_type
            case 1
                negViewEmo(1:20) = 1;
                posRegSuccess(21:40) = 1; 
                posRegFail(21:40) = 1; 
                neuView(41:60) = 1;
                negRegSuccess(61:80) = 1;
                negRegFail(61:80) = 1;
                posViewEmo(81:100) = 1;
            case 2
                posViewEmo(1:20) = 1;
                negRegSuccess(21:40) = 1;
                negRegFail(21:40) = 1;
                neuView(41:60) = 1;
                posRegSuccess(61:80) = 1;
                posRegFail(61:80) = 1;
                negViewEmo(81:100) = 1;
            case 3
                posRegSuccess(1:20) = 1;
                posRegFail(1:20) = 1;
                negViewEmo(21:40) = 1;
                neuView(41:60) = 1;
                posViewEmo(61:80) = 1;
                negRegSuccess(81:100) = 1;
                negRegFail(81:100) = 1;
            case 4
                negRegSuccess(1:20) = 1;
                negRegFail(1:20) = 1;
                posViewEmo(21:40) = 1;
                neuView(41:60) = 1;
                negViewEmo(61:80) = 1;
                posRegSuccess(81:100) = 1;
                posRegFail(81:100) = 1;
            case 5
                posViewEmo(1:20) = 1;
                negRegSuccess(21:40) = 1;
                negRegFail(21:40) = 1;
                neuView(41:60) = 1;
                negViewEmo(61:80) = 1;
                posRegSuccess(81:100) = 1;
                posRegFail(81:100) = 1;
        end
        
        % read in ratings - ratings only work for non-missed trials
        for i = 1 : nrTrials
            
            if ~isnan(ss.responses.responses(i))
                % er = emotion rating of this picture: parametric modulator
                % contains ratings for every trial
                er(i) = ss.responses.responses(i);
                rts(i) = ss.timing.reaction_times(i);
            end
            
            
            if isnan(ss.responses.responses(i))
                % if rating was missed, fill in average rating for the block
                if i <= 20
                    er(i) = nanmean(ss.responses.responses(1:20));
                end
                if (i > 20) && (i <= 40)
                    er(i) = nanmean(ss.responses.responses(21:40));
                end
                if (i > 40) && (i <= 60)
                    er(i) = nanmean(ss.responses.responses(41:60));
                end
                if (i > 60) && (i <= 80)
                    er(i) = nanmean(ss.responses.responses(61:80));
                end
                if (i > 80) && (i <= 100)
                    er(i) = nanmean(ss.responses.responses(81:100));
                end
            end
            
        end %end trial loop
        
        % calculate emotion regulation success
        % default case: positive regulation success 
        % should be scanner < post rating
        er_reg_success = post_rating_reg - er; 
        
        % multiply the negative regulate blocks by -1 
        % to get negative regulation success: 
        % should be scanner > post rating
        switch imglist.rand_type
            case 1
                er_reg_success = cat(1, er_reg_success(1:60), ...
                    er_reg_success(61:80)*(-1), er_reg_success(81:end));
            case 2
                er_reg_success = cat(1, er_reg_success(1:20), ...
                    er_reg_success(21:40)*(-1), er_reg_success(41:end));
            case 3
                er_reg_success = cat(1, er_reg_success(1:80), ...
                    er_reg_success(81:100)*(-1));
            case 4
                er_reg_success = cat(1, er_reg_success(1:20)*(-1), ...
                    er_reg_success(21:end));
            case 5
                er_reg_success = cat(1, er_reg_success(1:20), ...
                    er_reg_success(21:40)*(-1), er_reg_success(41:end));
        end
        
        % classify emotion regulation failures and successes:
        for i = 1 : nrTrials
            
            % if this was a positive regulation trial and they failed to change
            % their feeling:
            if (er_reg_success(i)) <= 0 && (posRegSuccess(i) == 1)
                posRegSuccess(i) = 0;
            end
            
            % if this was a positive regulation trial and they managed to change
            % their feeling:
            if (er_reg_success(i) > 0) && (posRegFail(i) == 1)
                posRegFail(i) = 0;
            end
            
            % if this was a negative regulation trial and they failed to change
            % their feeling:
            if (er_reg_success(i) <= 0) && (negRegSuccess(i) == 1)
                negRegSuccess(i) = 0;
            end
            
            % if this was a negative regulation trial and they managed to change
            % their feeling:
            if (er_reg_success(i) > 0) && (negRegFail(i) == 1)
                negRegFail(i) = 0;
            end
            
            
        end
        
        
        
        %% set up regressors of interest
        
        % subtract the absolute start_time (trigger)
        % from all time stamps, so the time line of the onsets starts at 0
        ss.timing.trial_start_times = ss.timing.trial_start_times-ss.timing.start_time; %start time for each choice
        ss.timing.content_start_times = ss.timing.content_start_times-ss.timing.start_time;
        
        
        % onsets
        ons{1,session} = ss.timing.content_start_times(posViewEmo > 0); % image to view / regulate comes on
        ons{2,session} = ss.timing.content_start_times(posRegSuccess > 0);
        ons{3,session} = ss.timing.content_start_times(posRegFail > 0);
        ons{4,session} = ss.timing.content_start_times(negViewEmo > 0);
        ons{5,session} = ss.timing.content_start_times(negRegSuccess > 0);
        ons{6,session} = ss.timing.content_start_times(negRegFail > 0);
        ons{7,session} = ss.timing.content_start_times(neuView > 0);
        ons{8,session} = ss.timing.trial_start_times; % rating screen comes on
        
        % durations  - as vectors
        ons_duration{1,session} = (ss.timing.trial_start_times(posViewEmo > 0)-ss.timing.content_start_times(posViewEmo > 0))'; % around 7 sec: from start content onset to start rating screen
        ons_duration{2,session} = (ss.timing.trial_start_times(posRegSuccess > 0)-ss.timing.content_start_times(posRegSuccess > 0))';
        ons_duration{3,session} = (ss.timing.trial_start_times(posRegFail > 0)-ss.timing.content_start_times(posRegFail > 0))';
        ons_duration{4,session} = (ss.timing.trial_start_times(negViewEmo > 0)-ss.timing.content_start_times(negViewEmo > 0))'; 
        ons_duration{5,session} = (ss.timing.trial_start_times(negRegSuccess > 0)-ss.timing.content_start_times(negRegSuccess > 0))';
        ons_duration{6,session} = (ss.timing.trial_start_times(negRegFail > 0)-ss.timing.content_start_times(negRegFail > 0))';
        ons_duration{7,session} = (ss.timing.trial_start_times(neuView > 0)-ss.timing.content_start_times(neuView > 0))'; 
        ons_duration{8,session} = ss.timing.reaction_times'; % all reaction times for the ratings - might contain NaN's
        % replace the NaNs with 4 seconds, screen was on for the rating time
        Replace(ons_duration{8,1}, NaN,4);
        
        
        
        %% movement - find bad onsets and define separate regressor for 3 TRs before and 1 after
        
        clear mov_reg{session}; % careful! do NOT clear mov_reg as a whole or will remove needed data for session > 1
        
        rpfile = dir([H_folder '/rp_*.txt']);
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
        
        if length(xcsmovindex)>=1 %just do this if there is any excess movement
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
    
    % fill information for the parametric modulators (pmod) for each onset
    pmod(1).name = {''};  % no modulator 
    pmod(1).param = {[]}; % actual values of the parametric modulator
    pmod(1).poly = {1}; % polynomial expansions, keep on 1 unless you want also quadratic effects
    
    pmod(2).name = {''};  % no modulator 
    pmod(2).param = {[]};
    pmod(2).poly = {1};
    
    pmod(3).name = {''}; % no modulator 
    pmod(3).param = {[]};
    pmod(3).poly = {1};
    
    pmod(4).name = {''}; % no modulator 
    pmod(4).param = {[]};
    pmod(4).poly = {1};
    
    pmod(5).name = {''}; % no modulator 
    pmod(5).param = {[]};
    pmod(5).poly = {1};
    
    pmod(6).name = {''}; % no modulator 
    pmod(6).param = {[]};
    pmod(6).poly = {1};
    
    pmod(7).name = {''}; % no modulator 
    pmod(7).param = {[]};
    pmod(7).poly = {1};
    
    pmod(8).name = {''}; % no modulator 
    pmod(8).param = {[]};
    pmod(8).poly = {1};
    
    
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
            t = t+1;
        else
            % write out
            low_event_matrix{temp5,1} = sub;
            low_event_matrix{temp5,2} = names{p};
            low_event_matrix{temp5,3} = length(onsets{p});
            temp5 = temp5+1;
        end
    end;
    
    
    % general parameters / SPM defaults for first level model
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
    number_of_TRs = ndyns; % length(scans)
    matlabbatch{m}.spm.stats.fmri_spec.sess(1).scans = scans;
    
    % load physio
    matlabbatch{m}.spm.stats.fmri_spec.sess(1).multi_reg = {reg_file_name};

    
    %% Module 2 for the SPM batch: save and estimate the SPM
    
    m = m+1;
    
    spm_jobman('initcfg');
    save_path2=[save_path filesep sub];
    
    matlabbatch{m}.spm.stats.fmri_est.spmmat = {[save_path2 filesep 'SPM.mat']};
    matlabbatch{m}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{m}.spm.stats.fmri_est.method.Classical = 1; % use Restricted Maximum Likelihood (REML) method
    
    spm_jobman('run',matlabbatch)
    
    clear matlabbatch
    
end

%% First level contrasts (start over with a new batch, start again at module 1)

m = 1;

% contrasts
CON_NAMES = {    'RegSuccess-ViewEmoME'...
    };

% These two lists must have the same length. If you have 1 sample ttests, then use something like
% null for the negative conname (can be any name that will not match one of the regressors)
connames_pos = { 'RegSuccess*bf',   ...
    };
connames_neg = {'ViewEmo*bf',  ...
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
    matlabbatch{m}.spm.stats.con.consess{i}.tcon.sessrep = 'none';  %don't replicate over sessions (replsc) if using create_cons function b/c already done
end
matlabbatch{m}.spm.stats.con.delete = 1; % write out all contrasts again if do cononly


% execute job and create con.nii files
spm_jobman('run', matlabbatch)

save([save_path2 filesep 'con_names_' sub], 'CON_NAMES') % save images for re-labeling contrasts later if cons get dropped bc not populated
save model_cons CON_NAMES; % save overview of all cons in the model results folder

    function con = create_cons(connames_pos, connames_neg, SPM)  % adapted from Jan Glaescher's code (I think)
        
        temp = load(SPM);
        SPM = temp.SPM;
        
        conweight = 1;
        negweight = -1;
        
        con = zeros(length(connames_pos),length(SPM.xX.name));
        negcon = zeros(length(connames_pos),length(SPM.xX.name));
          
        
        for c = 1:length(connames_pos)-1
            
            for r = 1:length(SPM.xX.name)
                if regexp(SPM.xX.name{r}, regexptranslate('wildcard',connames_pos{c}))
                    con(c,r) = conweight;
                end
            end
            
            
            for r = 1:length(SPM.xX.name)
                if regexp(SPM.xX.name{r}, regexptranslate('wildcard',connames_neg{c}))
                    negcon(c,r) = negweight;
                end
            end
            
            % Normalize contrast weights to account for inequalities across sessions
            con(c,:) = con(c,:) ./ sum(con(c,:));
            negcon(c,:) = -1*(negcon(c,:) ./ min(-1,sum(negcon(c,:)))); %Multiplying by -1 preserves the sign, min(-1,sum) prevents division by 0
            con(c,:) = con(c,:) + negcon(c,:);
            
        end % end normal con weights loop
        
    end % end con function

end