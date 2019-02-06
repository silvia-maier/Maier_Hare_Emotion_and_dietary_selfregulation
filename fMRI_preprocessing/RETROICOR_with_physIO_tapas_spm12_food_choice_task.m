% Create batch to run RETROICOR with the TAPAS physIO toolbox 
% integrated in SPM12 (6225)

% Author: Silvia Maier, University of Zurich, Zurich Center for Neuroeconomics 
% Contact: silvia.maier@econ.uzh.ch
% based on code contributions by Micah Edelson

% This script requires that the data were already preprocessed with 
% preprocessing_food_choice.m so that we here can retrieve the motion 
% realignment parameters

clear all; close all;

% select participants to preprocess
subs = {'SNS_ESC1_1003', 'SNS_ESC1_1004', ...
     'SNS_ESC1_1006', 'SNS_ESC1_1007', 'SNS_ESC1_1008', 'SNS_ESC1_1009', 'SNS_ESC1_1012', ...
     'SNS_ESC1_1014', 'SNS_ESC1_1015', 'SNS_ESC1_1016', ...
     'SNS_ESC1_1018', 'SNS_ESC1_1019', 'SNS_ESC1_1020', ...
     'SNS_ESC1_1021', 'SNS_ESC1_1023', 'SNS_ESC1_1026', 'SNS_ESC1_1027', ...
     'SNS_ESC1_1030', 'SNS_ESC1_1032', 'SNS_ESC1_1033', 'SNS_ESC1_1034', ...
     'SNS_ESC1_1036', 'SNS_ESC1_1037', 'SNS_ESC1_1038', 'SNS_ESC1_1039', ...
     'SNS_ESC1_1040', 'SNS_ESC1_1043', 'SNS_ESC1_1046', 'SNS_ESC1_1047', ...
     'SNS_ESC1_1048', 'SNS_ESC1_1049', 'SNS_ESC1_1050', 'SNS_ESC1_1052', ...
     'SNS_ESC1_1054', 'SNS_ESC1_1060', 'SNS_ESC1_1062', 'SNS_ESC1_1066', ...
     'SNS_ESC1_1067', 'SNS_ESC1_1068', 'SNS_ESC1_1071', 'SNS_ESC1_1072'
    };

% run the script directly without opening the SPM batch editor
run_script = 1;
open_batch = 0;

for i = 1:length(subs)
    
        % paths to raw data
        H_folder = ['/Users/username/mrdata/ESC/hrdata/' subs{i} filesep 'fc/']; % folder in which the physio files from the scanner are stored
        D_folder = ['/Users/username/mrdata/ESC/subs/' subs{i} filesep 'functional/fc/']; % folder with the functional scans for this task, fc = food choice
        OP_folder =['/Users/username/mrdata/ESC/hrdata/' subs{i} filesep 'fc/physio_output/']; % where the output is saved
        mkdir(OP_folder)
        
        % path to files with ECG and breathing belt information
        Physio_file_ID = dir([H_folder '*SCANPHYSLOG*']); % name of the physio files acquired with the Philips Achieva scanner in the SNS Lab
        Physio_path_and_file_ID = fullfile(H_folder, Physio_file_ID.name); 
        
        % create matlabbatch
        matlabbatch{1}.spm.tools.physio.save_dir = {OP_folder};
        matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Philips';
        matlabbatch{1}.spm.tools.physio.log_files.cardiac = {Physio_path_and_file_ID}; % location of the log file for cardiac information
        matlabbatch{1}.spm.tools.physio.log_files.respiration = {Physio_path_and_file_ID}; % same path for respiration (breathing belt) data
        matlabbatch{1}.spm.tools.physio.log_files.scan_timing = [];
        matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = 0.002; % 2 millisecond time resolution for the acquisition of ECG and breathing data
        matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
        % below specify the alignment of the scanner data to the log file - 
        % "last" means that the physIO toolbox will start looking from 
        % the last volume and count back the indicated number of volumes 
        % that you acquired; this is our best option because the physIO recording 
        % starts already during scanner preparation, before the actual run
        matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'last'; 
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = 37; % how many slices were acquired
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = []; % only for triggered (gated) sequences
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = 2.344; % repetition time
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 5; % number of dummies - how many discarded volumes were recorded in the beginning to allow the field to stabilize
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = 354; % number of volumes that were acquired within the run
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = 18; % reference slice - put the same as in the preprocessing (refslice = floor(nr slices/2))
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = []; % time between slices; if empty set as defult to TR/N slices
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = []; % count preparation pulses before 1 dummy - only fill if you use "first scan"
        matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.grad_direction = 'y'; % use this gradient direction to detect heartbeat signal
        matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.zero = 0.4;
        matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.slice = 0.45;
        matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.vol = []; % leave [] if unused; set value >= slice, if volume start gradients are higher than slice gradients
        matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.vol_spacing = []; % leave [] if unused
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'ECG';
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]); % allows to manually correct missing data (<20 pulses missing)
        matlabbatch{1}.spm.tools.physio.model.type = 'RETROICOR';
        matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
        matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
        matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
        matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
        % below commented out because we modeled movement separately
%         matlabbatch{1}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {[D_folder,'\rp_sn_*.txt']};
%         matlabbatch{1}.spm.tools.physio.model.movement.yes.outlier_translation_mm = Inf;
%         matlabbatch{1}.spm.tools.physio.model.movement.yes.outlier_rotation_deg = Inf;
        matlabbatch{1}.spm.tools.physio.model.input_other_multiple_regressors = {[D_folder,'\rp_sn_*.txt']};
        matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = ['fc_multiple_regressors.mat'];% output file
        matlabbatch{1}.spm.tools.physio.verbose.level = 2; % verbosity level of figures (can enter 1 (no output) until 3 for intense debugging)
        matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = 'PhysIO_output_level2.fig'; % output figure names
        matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;
        
        %%%run batch
        if run_script==1
            spm_jobman('run',matlabbatch)
            clear matlabbatch
        elseif open_batch==1
            spm_jobman('interactive',matlabbatch)
        end
        
    close all
end