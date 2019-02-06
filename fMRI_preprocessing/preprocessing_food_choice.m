% Preprocessing batch script for the food choice task of Maier &
% Hare "Emotion and dietary self-regulation"

% Author: Silvia Maier, University of Zurich, Zurich Center for Neuroeconomics
% Contact: silvia.maier@econ.uzh.ch
% based on code contributions by Aidan Makwana and Micah Edelson

%% Choose preprocessing steps to execute

% The following preprocessing steps were chosen for this project:
% (0 = no; 1 = yes)

toRun.realign = 0; % to only run realign: set to 1 if no fieldmap was acquired
% if fieldmap was acquired: (yes, was the case in this project)
toRun.fieldmap = 1; % run fieldmap preprocessing step
toRun.applyVDM = 0; % do not apply VDM
toRun.realignNunwarp = 1; % use realign & unwarp instead
toRun.sliceTiming = 1; % use slice timing
% warping steps:
% 1: if you want to use EPI template
toRun.EPIcoregister = 0;
toRun.EPInormalise = 0;
% 2: if you want to use the T1 template
toRun.T1coregister = 1;
toRun.T1normalise = 0;
% 3: if you want to use newsegment
toRun.NewSegment = 1;
toRun.deformations = 1;
toRun.smooth = 1;


%% Start spm_jobman
spm_jobman('initcfg');

%% Set paths and parameters

% set up directory names
data_dir_prefix = '/Users/username/mrdata/ESC';

folder.base = [data_dir_prefix]; % study level folder
folder.storage = [data_dir_prefix]; % where output data should be stored
folder.subject ='subs'; % folder with participant data
folder.parFiles ='parFiles'; % .par folder path from subject folder
folder.fieldmap ='fieldmap'; % folder path from participant level to your B0 maps
folder.anatomical ='anatomical'; % folder with structural T1 scan
folder.functional ='functional/fc/'; % folder with the functional images for this task (fc = food choice)
folder.spm12 ='/Applications/spm12/'; % where SPM lives
folder.scripts ='/Users/username/mrscripts/esc/preprocessing'; % folder where matlabbatches for preprocessing should be saved

ref.EPI = [folder.spm12 filesep 'templates/EPI.nii']; % for standalone SPM. this version won't work for normal SPM.
ref.T1 = [folder.spm12 filesep 'templates/T1.nii,1'];

% In this project, all short and long TEs for the B0 map
% for all participants are the same, set as below:
if toRun.fieldmap == 1 || toRun.applyVDM == 1 || toRun.realignNunwarp == 1
    searchString.B0 = '*b0*'; % assumes b0 order is in session order
    suffix.shtPhs = '_ec1_typ3.nii';
    suffix.shtMag = '_ec1_typ0.nii';
    suffix.lngPhs = '_ec2_typ3.nii';
    suffix.lngMag = '_ec2_typ0.nii';
end
shortTE = 4.299; % short echo time for the B0 map
longTE = 7.4; % long echo time for the B0 map


% Define a regexp search string to find participant data in the "subs" folder using 'searchString.subject'
searchString.subject = 'SNS_ESC1_[0-9][0-9][0-9][0-9]*'; % search string to get participants
% Define a dir search string to find all functional files in the functional folder using 'searchString.EPI'
searchString.EPI = '*run*'; % should be uniquely distinguishable from other file types
% Define a dir search string to find the T1 image in the subs folder using 'searchString.T1'
searchString.T1 = '*t1w3*'; % should be uniquely distinguishable from other file types

folder.subjects = dir([folder.storage filesep folder.subject]);
folder.subjects([folder.subjects.isdir] == 0) = [];
folder.subjects(cellfun(@isempty,regexp({folder.subjects.name},searchString.subject))) = [];
% subjects = {folder.subjects.name}; % or can enter manually as below
subjects = {'SNS_ESC1_1003', 'SNS_ESC1_1004', ...
    'SNS_ESC1_1006', 'SNS_ESC1_1007', 'SNS_ESC1_1008', 'SNS_ESC1_1009', 'SNS_ESC1_1012', ...
    'SNS_ESC1_1013', 'SNS_ESC1_1014', 'SNS_ESC1_1015', 'SNS_ESC1_1016', ...
    'SNS_ESC1_1017', 'SNS_ESC1_1018', 'SNS_ESC1_1019', 'SNS_ESC1_1020', ...
    'SNS_ESC1_1021', 'SNS_ESC1_1023', 'SNS_ESC1_1026', 'SNS_ESC1_1027', ...
    'SNS_ESC1_1030', 'SNS_ESC1_1032', 'SNS_ESC1_1033', 'SNS_ESC1_1034', ...
    'SNS_ESC1_1036', 'SNS_ESC1_1037', 'SNS_ESC1_1038', ...
    'SNS_ESC1_1039', ...
    'SNS_ESC1_1040', 'SNS_ESC1_1043', 'SNS_ESC1_1046', 'SNS_ESC1_1047', ...
    'SNS_ESC1_1048', 'SNS_ESC1_1049', 'SNS_ESC1_1050', 'SNS_ESC1_1052', ...
    'SNS_ESC1_1054', 'SNS_ESC1_1060', 'SNS_ESC1_1062', 'SNS_ESC1_1066', ...
    'SNS_ESC1_1067', 'SNS_ESC1_1068', 'SNS_ESC1_1071', 'SNS_ESC1_1072'
    };

% set prefixes to be used for output of each preprocessing step:
prefix.preProcVersion = '';
prefix.unwarp = 'u';
prefix.normalise = 'w';
prefix.smooth = 's';
prefix.T1coreg = 'cor';
prefix.realign = 'r';
prefix.sliceTiming = 'a';
prefix.currPrefix = ''; % if you wanted to start the preprocessing halfway through, you would need to put in the relevant prefixes i.e. 'upp1_'
% prefix.unwarp = ''; % uncomment if you don't want to use fieldmaps

%% setup study specifics

% set up runs
runs = 1; % all conditions of the paradigm were acquired in one single run

% set parameters for preprocessing
blip_d = -1; %blip direction: negative = -1; positive = 1
pullback = 1; %for deformations of type "pullback" enter 1; for "push forward" enter 0
defor_f_anat = 1; %run deformations from anatomical or functional segmentation

% slice timing: enter number of slices
options.slice_timing.Nslices = 37;

% EPI sequence
options.EPI.TR = 2.344; %repetition time (TR) in seconds
 % voxel resolution in millimeters (here: 2.5 x 2.5 x 3); 
 % for the slice thickness (z-direction), also add the slice gap (here: 0.6mm)
 % (open .par file with text editor to see header information 
 % in the matrix columns under "thick" and "gap")
options.EPI.nat_res = [2.5 2.5 3.6]; 
options.EPIreadoutTime = 39.0685; %look up the parameters in the header 
% information and calculate as below:

% Formula according to the FSL mailing list and using the "Philips Achieva
% constant" of 434.215 with our SNS Lab scanner according to Roger Luechinger:

% ( (1000 * Water Fat Shift[pixels]) / (434.215 * EPI factor + 1) ) * (EPI factor - 1)

% FSL mailing list: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;162ab1a3.1308


%% Build preprocessing batch

k = 0;

for i = 1:length(subjects)
    
    % Read in data for preprocessing
    
    fieldmapPath = [folder.base filesep folder.subject filesep subjects{i} filesep folder.fieldmap];
    anatPath = [folder.base filesep folder.subject filesep subjects{i} filesep folder.anatomical];
    funcPath = [folder.base filesep folder.subject filesep subjects{i} filesep folder.functional];
    parPath = [folder.base filesep folder.subject filesep subjects{i} filesep folder.parFiles];
    
    % retrieve T1 scan
    try
        T1anat = dir([anatPath filesep subRef searchString.T1 'nii']);
    end
    
    % pull out all functional scans
    subRef = 'sn_'; % all raw functional images start with this prefix
    subsearchString.EPI = [subRef searchString.EPI 'nii'];
    scannames = dir([funcPath filesep subsearchString.EPI]);
    
    try
        for r = 1:length(runs)
            scans{r} = {scannames(r).name};
        end
    end
    
    for r = 1:length(runs)
        if strcmp(prefix.currPrefix,'')
            currPrefix=prefix.preProcVersion;
            for n = 1:length(scans{r});
                try
                    copyfile([funcPath filesep scans{r}{n}],[funcPath filesep currPrefix scans{r}{n}]);
                end;
            end
            try
                copyfile([anatPath filesep T1anat.name],[anatPath filesep prefix.preProcVersion T1anat.name]);
            end
        else
            currPrefix=prefix.currPrefix;
        end
        
        for n = 1:length(scans{r});
            scanFiles{r}{n} = cellstr(spm_select('ExtList',funcPath, ['^' scans{r}{n}],inf));
        end
    end
    
    %% FIELDMAP (PHASE AND MAG)
    
    if toRun.fieldmap == 1
        B0maps = dir([fieldmapPath filesep subRef searchString.B0 'nii']);
        
        for n = 1:length(scans)
            
            B0mapInd = 1;  % use run 1 phase and mag, we just have 1 run in this project
            B0map = B0maps(B0mapInd).name(1:end-13);
            shortPhs = [fieldmapPath filesep B0map suffix.shtPhs ',1'];
            shortMag = [fieldmapPath filesep B0map suffix.shtMag ',1'];
            longPhs = [fieldmapPath filesep B0map suffix.lngPhs ',1'];
            longMag = [fieldmapPath filesep B0map suffix.lngMag ',1'];
            k = k+1;
            
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.shortphase = {shortPhs};
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.shortmag = {shortMag};
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.longphase = {longPhs};
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.longmag = {longMag};
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.et = [shortTE longTE]; % short TE and long TE of B0
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.maskbrain = 1;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.blipdir = blip_d;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.tert = options.EPIreadoutTime;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.epifm = 0;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.ajm = 0;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.method = 'Mark3D';
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.pad = 0;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.ws = 1;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.template = {'/Applications/spm12/toolbox/FieldMap/T1.nii'};
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.nerode = 2;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.thresh = 0.5;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.reg = 0.02;
            for r = 1:length(runs)
                matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.session(r).epi = {[funcPath filesep currPrefix scanFiles{r}{:}{1}]};
            end
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.matchvdm = 1;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.sessname = 'session';
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.writeunwarped = 1;
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.anat = {''};
            matlabbatch{k}.spm.tools.fieldmap.phasemag.subj.matchanat = 1;
        end
    end
    
    %% REALIGN & UNWARP
    
    if toRun.realignNunwarp == 1
        k = k+1;
        
        for r = 1:length(runs);
            matlabbatch{k}.spm.spatial.realignunwarp.data(r).scans = strcat([funcPath filesep currPrefix], scanFiles{r}{:});
            matlabbatch{k}.spm.spatial.realignunwarp.data(r).pmscan = ''; % if left blank, spm will estimate b0
        end
        
        %%%put in parameters
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.rtm = 0; % register to mean (1) or to first (0) image
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.einterp = 7; % level of interpolation; have same as below in rinterp!
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
        matlabbatch{k}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{k}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{k}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{k}.spm.spatial.realignunwarp.uwroptions.rinterp = 7; % level of interpolation
        matlabbatch{k}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
        matlabbatch{k}.spm.spatial.realignunwarp.uwroptions.mask = 1;
        matlabbatch{k}.spm.spatial.realignunwarp.uwroptions.prefix = prefix.unwarp;
        currPrefix=[prefix.unwarp currPrefix];
    end
    
    
    %% SLICE TIMING
    
    if toRun.sliceTiming == 1
        k = k+1;
        for r = 1:length(runs)
            matlabbatch{k}.spm.temporal.st.scans{r} = strcat([funcPath filesep currPrefix],scanFiles{r}{:});
        end
        matlabbatch{k}.spm.temporal.st.nslices =  options.slice_timing.Nslices;
        matlabbatch{k}.spm.temporal.st.tr = options.EPI.TR;
        matlabbatch{k}.spm.temporal.st.ta = options.EPI.TR-options.EPI.TR/options.slice_timing.Nslices;
        matlabbatch{k}.spm.temporal.st.so = 1:options.slice_timing.Nslices; % ascending slice acquisition
        matlabbatch{k}.spm.temporal.st.refslice = floor(options.slice_timing.Nslices/2);  % round for even number of slices, floor or ceiling for uneven number % /3 if OFC centered
        matlabbatch{k}.spm.temporal.st.prefix = prefix.sliceTiming;
        currPrefix = [prefix.sliceTiming currPrefix];
    end
    
    
    %% EPI COREGISTER
    
    if toRun.EPIcoregister == 1
        k = k+1;
        matlabbatch{k}.spm.spatial.coreg.estimate.ref = {ref.EPI};
        matlabbatch{k}.spm.spatial.coreg.estimate.source = {[funcPath filesep 'mean' prefix.preProcVersion, scanFiles{1}{1}{1}]};
        matlabbatch{k}.spm.spatial.coreg.estimate.other = strcat([funcPath filesep currPrefix],vertcat(scanFiles{1}{:}));
        matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    end
    
    %% T1 COREGISTER
    
    if toRun.T1coregister == 1
        k = k+1;
        matlabbatch{k}.spm.spatial.coreg.estwrite.ref(1) = {[funcPath filesep 'meanu' prefix.preProcVersion, scanFiles{1}{1}{1}(1:end-2)]};
        matlabbatch{k}.spm.spatial.coreg.estwrite.source = strcat(anatPath, filesep, prefix.preProcVersion, {T1anat.name}, ',1');
        matlabbatch{k}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi'; % normalised mutual information
        matlabbatch{k}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{k}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; % spm default
        matlabbatch{k}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{k}.spm.spatial.coreg.estwrite.roptions.interp = 4; % get SPM default:  @(val)spm_get_defaults('coreg.write.interp',val{:})
        matlabbatch{k}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{k}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{k}.spm.spatial.coreg.estwrite.roptions.prefix = prefix.T1coreg;
    end
    
    %% EPI NORMALISE
    
    if toRun.EPInormalise == 1
        for n = 1:length(scans{r})
            k = k+1;
            matlabbatch{k}.spm.spatial.normalise.estwrite.subj.source = {[funcPath filesep 'mean' prefix.preProcVersion, scanFiles{1}{1}{1}]};
            matlabbatch{k}.spm.spatial.normalise.estwrite.subj.wtsrc = {};
            matlabbatch{k}.spm.spatial.normalise.estwrite.subj.resample = strcat([funcPath filesep currPrefix], vertcat(scanFiles{1}{:}));
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.template = {ref.EPI};
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.weight = {''};
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
            matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
            matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
            matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.bb =[-78 -112 -50; 78 76 85];
            matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.vox = options.EPI.nat_res;
            matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.interp = 4;
            matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.prefix = prefix.normalise;
        end
        currPrefix = [prefix.normalise currPrefix];
    end
    
    
    %% T1 NORMALISE
    
    if toRun.T1normalise == 1
        k = k+1;
        matlabbatch{k}.spm.spatial.normalise.estwrite.subj.source = strcat(anatPath, filesep, prefix.T1coreg, prefix.preProcVersion, {T1anat.name}, ',1');
        matlabbatch{k}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
        matlabbatch{k}.spm.spatial.normalise.estwrite.subj.resample = strcat([funcPath filesep currPrefix], vertcat(scanFiles{1}{:}));
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.template = {ref.T1};
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.weight = '';
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
        matlabbatch{k}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
        matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
        matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50; 78 76 85];
        matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.vox = options.EPI.nat_res;
        matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.interp = 1;
        matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{k}.spm.spatial.normalise.estwrite.roptions.prefix = prefix.normalise;
        currPrefix = [prefix.normalise currPrefix];
    end
    
    
    %% NEW SEGMENT
    
    % currently set to write forward + backward deformations for EPI and T1
    % and grey, white, csf, bone, soft tissue and air tissue classifications
    % (c1-6) for T1 only
    
    if toRun.NewSegment == 1
        
        % segment anatomical
        k = k+1;
        matlabbatch{k}.spm.spatial.preproc.channel.vols = strcat(anatPath, filesep, prefix.T1coreg, prefix.preProcVersion, {T1anat.name}, ',1');
        matlabbatch{k}.spm.spatial.preproc.channel.biasreg = 0.001; % light regularization, SPM default
        matlabbatch{k}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{k}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(1).tpm = {'/Applications/spm12/tpm/TPM.nii,1'};
        matlabbatch{k}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(1).native = [1 1]; % write out tissue class images (c*) native and/or dartel
        matlabbatch{k}.spm.spatial.preproc.tissue(1).warped = [1 1]; % write out (mwc*) and/or (wc*) files)
        matlabbatch{k}.spm.spatial.preproc.tissue(2).tpm = {'/Applications/spm12/tpm/TPM.nii,2'};
        matlabbatch{k}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(2).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(2).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(3).tpm = {'/Applications/spm12/tpm/TPM.nii,3'};
        matlabbatch{k}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(3).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(3).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(4).tpm = {'/Applications/spm12/tpm/TPM.nii,4'};
        matlabbatch{k}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{k}.spm.spatial.preproc.tissue(4).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(4).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(5).tpm = {'/Applications/spm12/tpm/TPM.nii,5'};
        matlabbatch{k}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{k}.spm.spatial.preproc.tissue(5).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(5).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(6).tpm = {'/Applications/spm12/tpm/TPM.nii,6'};
        matlabbatch{k}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(6).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(6).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.warp.mrf = 2; % markov random field cleanup
        matlabbatch{k}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{k}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]; % SPM default
        matlabbatch{k}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{k}.spm.spatial.preproc.warp.fwhm = 0; % smoothness 0 for MRI
        matlabbatch{k}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{k}.spm.spatial.preproc.warp.write = [1 1];
        
        % segment functional
        k = k+1;
        matlabbatch{k}.spm.spatial.preproc.channel.vols(1) = {[funcPath filesep 'meanu' prefix.preProcVersion, scanFiles{1}{1}{1}(1:end-2)]};
        matlabbatch{k}.spm.spatial.preproc.channel.biasreg = 0.001; % light regularization, SPM default
        matlabbatch{k}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{k}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(1).tpm = {'/Applications/spm12/tpm/TPM.nii,1'};
        matlabbatch{k}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(1).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(1).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(2).tpm = {'/Applications/spm12/tpm/TPM.nii,2'};
        matlabbatch{k}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(2).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(2).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(3).tpm = {'/Applications/spm12/tpm/TPM.nii,3'};
        matlabbatch{k}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(3).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(3).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(4).tpm = {'/Applications/spm12/tpm/TPM.nii,4'};
        matlabbatch{k}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{k}.spm.spatial.preproc.tissue(4).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(4).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(5).tpm = {'/Applications/spm12/tpm/TPM.nii,5'};
        matlabbatch{k}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{k}.spm.spatial.preproc.tissue(5).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(5).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(6).tpm = {'/Applications/spm12/tpm/TPM.nii,6'};
        matlabbatch{k}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(6).native = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(6).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.warp.mrf = 2;
        matlabbatch{k}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{k}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{k}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{k}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{k}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{k}.spm.spatial.preproc.warp.write = [1 1];
        
    end
    
    
    %% DEFORMATIONS
    
    if toRun.deformations == 1
        % deformations anatomical
        k = k+1;
        if pullback == 1 % pullback deform method
            matlabbatch{k}.spm.util.defs.comp{1}.def = strcat(anatPath, filesep,'y_',prefix.T1coreg, prefix.preProcVersion, {T1anat.name});
            matlabbatch{k}.spm.util.defs.out{1}.pull.fnames = strcat(anatPath, filesep,prefix.T1coreg, prefix.preProcVersion, {T1anat.name});
            matlabbatch{k}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
            matlabbatch{k}.spm.util.defs.out{1}.pull.interp = 7;
            matlabbatch{k}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{k}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            
        else  % pushforward
            matlabbatch{k}.spm.util.defs.comp{1}.def(1) = strcat(anatPath, filesep,'iy_',prefix.T1coreg, prefix.preProcVersion, {T1anat.name});
            matlabbatch{k}.spm.util.defs.out{1}.pull.fnames = strcat(anatPath, filesep,prefix.T1coreg, prefix.preProcVersion, {T1anat.name});
            matlabbatch{k}.spm.util.defs.out{1}.push.weight = {''};
            matlabbatch{k}.spm.util.defs.out{1}.push.savedir.savesrc = 1;
            matlabbatch{k}.spm.util.defs.out{1}.push.fov.bbvox.bb = [inf inf inf
                inf inf inf];
            matlabbatch{k}.spm.util.defs.out{1}.push.fov.bbvox.vox = [2.5 2.5 3.6];
            matlabbatch{k}.spm.util.defs.out{1}.push.preserve = 0;
            matlabbatch{k}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        end
        
        % deformations functional (8)
        
        k = k+1;
        if pullback == 1 % pullback deform method
            if defor_f_anat == 1
                matlabbatch{k}.spm.util.defs.comp{1}.def = strcat(anatPath, filesep,'y_',prefix.T1coreg, prefix.preProcVersion, {T1anat.name});
            else % change below if defor_f_anat ~= 1
                matlabbatch{k}.spm.util.defs.comp{1}.def = {[funcPath filesep 'y_meanu' prefix.preProcVersion, scanFiles{1}{1}{1}(1:end-2)]};
            end
            matlabbatch{k}.spm.util.defs.comp{2}.idbbvox.vox = [2.5 2.5 3.6];
            matlabbatch{k}.spm.util.defs.comp{2}.idbbvox.bb = [inf inf inf
                inf inf inf];
            
            for r = 1:length(runs)
                matlabbatch{1,k}.spm.util.defs.out{1, 1}.pull.fnames(r) = {[funcPath filesep currPrefix scanFiles{r}{1}{1}(1:end-2)]};
            end;
            
            matlabbatch{k}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
            matlabbatch{k}.spm.util.defs.out{1}.pull.interp = 7; % Anjali: 4
            matlabbatch{k}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{k}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            
        else % pushforward
            if defor_f_anat == 1
                matlabbatch{k}.spm.util.defs.comp{1}.def = strcat(anatPath, filesep,'iy_',prefix.T1coreg, prefix.preProcVersion, {T1anat.name});
            else
                matlabbatch{k}.spm.util.defs.comp{1}.def = {[funcPath filesep 'iy_meanu' prefix.preProcVersion, scanFiles{1}{1}{1}(1:end-2)]};
            end
            
            matlabbatch{k}.spm.util.defs.comp{2}.idbbvox.vox = [2.5 2.5 3.6];
            matlabbatch{k}.spm.util.defs.comp{2}.idbbvox.bb = [inf inf inf
                inf inf inf];
            
            for r = 1:length(runs)
                matlabbatch{1, k}.spm.util.defs.out{1, 1}.pull.fnames(r) = strcat([funcPath filesep currPrefix], scanFiles{r}{1}{1}(1:end-2));
            end;
            
            matlabbatch{k}.spm.util.defs.out{1}.push.weight = {''};
            matlabbatch{k}.spm.util.defs.out{1}.push.savedir.savesrc = 1;
            matlabbatch{k}.spm.util.defs.out{1}.push.fov.bbvox.bb = [inf inf inf
                inf inf inf];
            matlabbatch{k}.spm.util.defs.out{1}.push.fov.bbvox.vox = [2.5 2.5 3.6];
            matlabbatch{k}.spm.util.defs.out{1}.push.preserve = 0;
            matlabbatch{k}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        end
        
        % currPrefix=['w' currPrefix]; % apparently implicit 'w' prefix
    end
    
    %% SMOOTH
    
    if toRun.smooth == 1
        k = k+1;
        for r = 1:length(runs)
            matlabbatch{k}.spm.spatial.smooth.data(r) = {[funcPath filesep prefix.normalise currPrefix scanFiles{r}{1}{1}(1:end-2)]};
        end
        matlabbatch{k}.spm.spatial.smooth.fwhm = [4 4 4]; % apply as little smoothing as possible to later run FSL Randomise
        matlabbatch{k}.spm.spatial.smooth.dtype = 0;
        matlabbatch{k}.spm.spatial.smooth.im = 0;
        matlabbatch{k}.spm.spatial.smooth.prefix = prefix.smooth;
    end
    
    save([folder.scripts filesep prefix.preProcVersion subjects{i} '_fc_run_batch.mat'],'matlabbatch');
    
    % if you want to run the batch at this point already: uncomment below
    %
    %     %%%run batch
    %     if run_script == 1
    %         spm_jobman('run',matlabbatch)
    %         clear matlabbatch
    %     elseif open_batch == 1
    %         spm_jobman('interactive',matlabbatch)
    %     end
    %
    
    clearvars -except runs run toRun folder ref searchString prefix suffix subjects options shortTE longTE blip_d pullback defor_f_anat run_script open_batch i r data_dir_prefix data_dir_prefix
    k = 0;
    
end % loop over all participants