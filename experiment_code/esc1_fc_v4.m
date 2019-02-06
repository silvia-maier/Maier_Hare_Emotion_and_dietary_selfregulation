function [responses] = esc1_fc_v4(subject_id, scale_direction, useEyetracker)
% fMRI task that presents food choices - option for the participant is
% to eat the food on the screen or nothing (respond: "Yes / No to food on 
% the display")

% inputs: 
% subject_id as a string, eg '1001'
% scale_direction: 'l_to_r' for left to right or 'r_to_l' for right to left
% (negative -> positive answer modality)
% useEyetracker: logical (1 = eyetracking is recorded, 0 = not recorded)


if ~exist('useEyetracker', 'var')
    useEyetracker = true;
end

nrRatings = 100; % total number of trials
TR = 2.344;  % TR (repetition time) for acquisition of 1 volume
decision_time = 3; % seconds for making a choice
HideCursor;


%% OUTPUT FILE
results_file_name = sprintf('food_choices_%s', subject_id);

% Check to prevent overwriting previous data
A = exist([results_file_name '.mat'], 'file');
if A
    writeover = input('Filename already exists, do you want to overwrite? 1=yes, 0=no ');
else
    writeover = 2; % no messages necessary
end

switch writeover
    case 0
        subject_id = input('Enter proper subject ID as text string: ');
    case 1
        disp('Results file will be overwritten')
end


%% PREPARE KEYBOARD
KbName('UnifyKeyNames');
escape_key = KbName('ESCAPE');
choice = KbName('space');


switch scale_direction
    case 'l_to_r'
        choice_key_names = {'g', 'b'}; % g = yes (left), b = no (right)
    case 'r_to_l'
        choice_key_names = {'b', 'g'}; % b = yes (left), g = no (right)
end

trigger_name = {'t'};

choice_keys = KbName(choice_key_names);
trigger_key = KbName(trigger_name);
escape_key = KbName('ESCAPE');


switch scale_direction
    case 'l_to_r'
        decision_names = {'Ja', 'Nein'};
    case 'r_to_l'
        decision_names = {'Nein', 'Ja'};
end

switch scale_direction
    case 'l_to_r'
        decision_scale = '     JA                             NEIN';
    case 'r_to_l'
        decision_scale = '     NEIN                           JA';
end


%% PREPARE INSTRUCTIONS

instructions{1} = ['Sie werden gleich ENTSCHEIDEN, was Sie am Ende der Studie essen moechten.' ...
    '\n\n Sie antworten mittels der folgenden Skala: \n\n' decision_scale '\n\n ' ...
    'Wenn Sie das Nahrungsmittel auf dem Bildschirm essen moechten, druecken Sie "Ja". ' ...
    'Wenn Sie das Nahrungsmittel auf dem Bildschirm nicht moechten, druecken Sie "Nein" ' ...
    'zum aktuellen Bildschirm. \n\n Ihr ZIEL in der Entscheidungsrunde sollte es sein, sich ' ...
    'MOEGLICHST GESUND zu ernaehren. An dieses Ziel wird Sie zwischen den Entscheidungen ein ' ...
    'kleines Symbol erinnern. \n\n Bitte entscheiden Sie sich immer so, wie Sie es als richtig ' ...
    'empfinden. Wenn Sie soweit sind, druecken Sie eine Taste um zu beginnen.'];


%% LOAD PICTURES and NAMES
folder = [pwd '/food_images']; % folder with Swiss 180 food item set
folder_scr = [pwd '/food_images_scr']; % folder with all scrambled food images
folder_fix = [pwd '/fixation_item']; % folder with the health symbol image

% pull all images from the folders
filenames = dir([folder '/*.bmp']); 
filenames_scr = dir([folder_scr '/*.bmp']);

% load the design (in which order images are presented to this participant)
fc = load(['imagelist_' subject_id]);
image_order = fc.imgvect.blocks; % images to present for choice
jitter = fc.imgvect.jitter; % for timing the ITI


%% PREPARE DISPLAY
[windowHandle, rectangle] = Screen('OpenWindow', max(Screen('Screens')), [0.5 0.5 0.5]); 

% careful! display of choice screen has been optimized 
% for resolution 1280 x 768 on the screen of the scanner PC!


% cell for loading the fixation item
fixationImage = cell(1, 1);
fixationImage{1, 1} = double(imread([folder_fix filesep ...
    'fixation_item_gray.bmp' ]));


% cell for loading the images
ratingImages = cell(nrRatings, 1);

for currRatingImage = 1 : nrRatings
    ratingImages{currRatingImage, 1} = double(imread([folder filesep ...
        filenames(image_order(currRatingImage)).name ]));
end


% cell for loading the scrambled images
ratingImagesScr = cell(nrRatings, 1);

for currRatingImage = 1 : nrRatings
    ratingImagesScr{currRatingImage, 1} = double(imread([folder_scr filesep ...
        filenames_scr(image_order(currRatingImage)).name ]));
end


%% PREPARE OUTPUT VARIABLES

timing = struct; % timing variable
responses = NaN(nrRatings,1); % responses for food choices

timing.trial_start_times = NaN(1,length(1:nrRatings)); % for food choices
timing.decision_times = NaN(1,length(1:nrRatings)); % time the choice was logged
timing.reaction_times = NaN(1,length(1:nrRatings));  % reaction time for the choice
timing.scramble_start_times =  NaN(1,length(1:nrRatings)); % for pupil adaptation: record start time of the cue

timing.blank_start_times = NaN(nrRatings, 1); % record start times for ITI


%% PREPARE SCREEN SETTINGS

% determine center of the opened screen
[width, height] = Screen('WindowSize', max(Screen('Screens')));
centerX = floor(0.5 * width);
centerY = floor(0.5 * height);

% colors and text
txt_color = [200 200 200]; % off white
bg_color = GrayIndex(windowHandle, 0.5); % middle gray, 0 = black, 1 = white
wrapat_length = 60;
Screen('TextFont', windowHandle,'Arial');

txt_size_prompt = 18;
txt_size_task = 14;

% image height and width
% food images:
imgH = 298;
imgW = 398;

% frame the images for selection -
% position of these frames optimized for screen resolution 1280 x 768 
frame1 = [rectangle(3)*3/8 rectangle(4)*4.2/6 rectangle(3)*3.9/8 rectangle(4)*4.5/6]; % frame left
frame2 = [rectangle(3)*4/8 rectangle(4)*4.2/6 rectangle(3)*4.9/8 rectangle(4)*4.5/6]; % frame right


%% CREATE TEXTURES

% fixation item to display
fixationTexture = cell(1, 1);
fixationTexture{1, 1} = Screen('MakeTexture', ...
    windowHandle, fixationImage{1, 1});

% food pics to display
ratingTexture = cell(nrRatings, 1);
for currRatingImage = 1 : nrRatings
    ratingTexture{currRatingImage, 1} = Screen('MakeTexture', ...
        windowHandle, ratingImages{currRatingImage, 1});
end

% food pics to display scrambled
ratingTextureScr = cell(nrRatings, 1);
for currRatingImage = 1 : nrRatings
    ratingTextureScr{currRatingImage, 1} = Screen('MakeTexture', ...
        windowHandle, ratingImagesScr{currRatingImage, 1});
end


%% DISPLAY PICTURES FOR RATING
try
    
    %drawing the instructions
    Screen('TextSize', windowHandle, txt_size_prompt);
    for i=1:length(instructions) 
        
        Screen(windowHandle, 'FillRect', bg_color);
        DrawFormattedText(windowHandle, instructions{i}, 'center', 'center', txt_color, wrapat_length, [], [], 1.5);
        Screen(windowHandle, 'Flip');
        WaitSecs(.5);
        
        while 1
            [key_is_down, ~, key_code] = KbCheck;
            
            % flip the instructions with any choice key!
            if key_is_down && any(key_code(choice_keys))
                break;
            end
            
        end
        
    end
    
    Screen(windowHandle, 'FillRect', bg_color);
    
    DrawFormattedText(windowHandle, 'Warten auf Scanner', 'center', 'center', txt_color);
    Screen(windowHandle, 'Flip');
    
    FlushEvents();
    
    % wait for scanner trigger (expects a "t" input!)
    waiting_for_scanner = 1;
    
    while waiting_for_scanner
        
        [key_is_down, ~, key_code] = KbCheck;
        
        if key_is_down && any(key_code(trigger_key))
            
            waiting_for_scanner = 0;
            timing.start_time = GetSecs; % start of the scanner timing!
            
            if useEyetracker
                %send message to EyeLink recording file
                Eyelink('message', 'ExpStart');
                WaitSecs(0.01); % wait 10 milliseconds to send the message
            end
            
        end
        
    end
    
    Screen('TextSize', windowHandle, txt_size_task);
    Screen(windowHandle, 'FillRect', bg_color);
    % draw the fixation cross
    Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
    Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
    Screen(windowHandle, 'Flip');
        
    % initial fixation period of 2 * TR, 
    % subtract the 0.01 seconds wait time for the eyetracker message (line 269) 
    % (we always used the eyetracker while scanning)
    WaitSecs(2*TR - 0.01); 
    
    for currRatingImage = 1 : nrRatings % main draw loop
        
        % DISPLAY SCRAMBLED IMAGE for adapting pupil response
        Screen('DrawTexture', windowHandle, ratingTextureScr{currRatingImage, 1});
        % with fixation cross on top
        Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
        Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
        
        Screen('flip',windowHandle);
        loopTimer = tic();
        timing.scramble_start_times(currRatingImage) = GetSecs();
        WaitSecs(1); % 1 second to allow the pupil to adapt
        
        % CHOICE SCREEN
        Screen('TextSize', windowHandle, txt_size_task);
        Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1}); 
        
        % with fixation cross on top
        Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
        Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
        
        % display answer scale below
        DrawFormattedText(windowHandle, decision_scale, centerX-140, centerY+155, txt_color, wrapat_length, [], [], 1.5);
        
        Screen('flip',windowHandle);
        timing.trial_start_times(currRatingImage) = GetSecs;
        
        % read out answer: what do they want to eat?
        choice_made = false;
        flipscreen = false;
        while ~choice_made
            
            [key_is_down, ~, key_code] = KbCheck;
            
            if key_is_down && any(key_code(choice_keys))&& length(find(key_code))<2 % pressing only one key at a time
                down_key = find(key_code, 1);
                responses(currRatingImage) = find(down_key == choice_keys, 1);
                timing.decision_times(currRatingImage) = GetSecs;
                timing.reaction_times(currRatingImage) = timing.decision_times(currRatingImage) - timing.trial_start_times(currRatingImage);
                choice_made = true;
            end
      
            if key_is_down && key_code(escape_key)
                save(results_file_name, 'responses', 'timing', 'jitter');
                break
            end
            
            if (GetSecs - timing.trial_start_times(currRatingImage)) >= decision_time % 3 seconds decision time
                break
            end
            
        end % while loop choice
        
        % answer feedback
        if strcmp(scale_direction, 'l_to_r')
            
            if responses(currRatingImage) == 1
                Screen('TextSize', windowHandle, txt_size_task);
                % draw the answer feedback as unobtrusive as possible to minimally disturb the pupil
                Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1}); 
                % with fixation cross on top
                Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
                Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
                % display answer scale below
                DrawFormattedText(windowHandle, decision_scale, centerX-140, centerY+155, txt_color, wrapat_length, [], [], 1.5);
                
                % frame their response
                Screen('FrameRect', windowHandle, txt_color, frame1, 1);
                Screen(windowHandle, 'Flip');
                WaitSecs(0.1);
            end
            
            if responses(currRatingImage) == 2
                Screen('TextSize', windowHandle, txt_size_task);
                % draw the answer feedback as unobtrusive as possible to minimally disturb the pupil
                Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1}); 
                % with fixation cross on top
                Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
                Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
                % display answer scale below
                DrawFormattedText(windowHandle, decision_scale, centerX-140, centerY+155, txt_color, wrapat_length, [], [], 1.5);
                
                % frame their response
                Screen('FrameRect', windowHandle, txt_color, frame2, 1);
                Screen(windowHandle, 'Flip');
                WaitSecs(0.1);
            end
            
        end
        
        if strcmp(scale_direction, 'r_to_l')
            
            if responses(currRatingImage) == 2
                Screen('TextSize', windowHandle, txt_size_task);
                Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1}); 
                % with fixation cross on top
                Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
                Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
                % display answer scale below
                DrawFormattedText(windowHandle, decision_scale, centerX-140, centerY+155, txt_color, wrapat_length, [], [], 1.5);
                
                % frame their response
                Screen('FrameRect', windowHandle, txt_color, frame1, 1);
                Screen(windowHandle, 'Flip');
                WaitSecs(0.1);
                
            end
            
            if responses(currRatingImage) == 1
                Screen('TextSize', windowHandle, txt_size_task);
                Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1});
                % with fixation cross on top
                Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
                Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
                % display answer scale below
                DrawFormattedText(windowHandle, decision_scale, centerX-140, centerY+155, txt_color, wrapat_length, [], [], 1.5);
                
                % frame their response
                Screen('FrameRect', windowHandle, txt_color, frame2, 1);
                Screen(windowHandle, 'Flip');
                WaitSecs(0.1);
            end
            
        end
        
        
        % fixation item - jittered ITI
        Screen('TextSize', windowHandle, txt_size_task);
        Screen(windowHandle, 'FillRect', bg_color);
        Screen('DrawTexture', windowHandle, fixationTexture{1,1},[], [rectangle(4)*4.25/9+(rectangle(3)-rectangle(4))/2  rectangle(4)*4.25/9 rectangle(4)*4.75/9+(rectangle(3)-rectangle(4))/2 rectangle(4)*4.75/9]);
        Screen('flip', windowHandle);
        timing.blank_start_times(currRatingImage) = GetSecs;
        
        save(results_file_name, 'responses',  'timing', 'jitter'); % save here to update results on every loop
        
        
        if useEyetracker
            %send message to EyeLink recording file
            Eyelink('message', sprintf('trial: %d', currRatingImage));
            Eyelink('command', sprintf('record_status_message ''End trial %d''', currRatingImage));
        end
        
        % cue present time + pic present time + jitter
        while toc(loopTimer) < 1 + 3 + jitter(currRatingImage), end
        timing.loopTimer(currRatingImage) = toc(loopTimer);
        
    end % end of main draw loop
    
    if useEyetracker
        % send message to EyeLink recording file
        Eyelink('message', 'ExpEnd');
        Eyelink('command', 'record_status_message ''ExpEnd''');
    end
    
    % fixation at the end of the run
    Screen('TextSize', windowHandle, txt_size_task);
    Screen(windowHandle, 'FillRect', bg_color);
    Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
    Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);    Screen('flip', windowHandle);
    WaitSecs(10*TR);

    save(results_file_name, 'responses',  'timing');
    
catch
    
    psychrethrow(psychlasterror);
    
end % end of try-catch


%% CLOSE DISPLAY
Screen('CloseAll');

end