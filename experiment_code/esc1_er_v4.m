function [responses] = esc1_er_v4(subject_id, scale_direction, useEyetracker)
% fMRI emotion reappraisal task that presents IAPS pictures

% inputs: 
% subject_id as a string, eg '1001'
% scale_direction: 'l_to_r' for left to right or 'r_to_l' for right to left
% (negative -> positive directionality)
% useEyetracker: logical (1 = eyetracking is recorded, 0 = not recorded)

% cue is presented before block, on top of the image just a small
% reminder to not divert gaze away from center of the image


if ~exist('useEyetracker', 'var')
    useEyetracker = true;
end

nrRatings = 100; % total number of trials
TR = 2.344; % TR (repetition time) for acquisition of 1 volume
display_time = 7; % seconds for viewing / reappraising
rating_time = 4; % seconds for making an emotion rating
HideCursor;
break_time = 15; % in seconds; time for pause between blocks


%% OUTPUT FILE
results_file_name = sprintf('er_scan_%s', subject_id);

% Check to prevent overwriting previous data
A = exist([results_file_name '.mat'], 'file');
if A
    writeover = input('Filename already exists, do you want to overwrite? 1 = yes, 0 = no ');
else
    writeover = 2; % no messages necessary
end

switch writeover
    case 0
        subject_id = input('Enter proper subject ID as text string: ');
    case 1
        disp('Results file will be overwritten')
end


%% PREPARE INSTRUCTIONS

instructions{1} = ['Sie werden in dieser Aufgabe eine Reihe von BILDERN gezeigt bekommen.' ...
    'Nach einem Block mit 20 Bildern gibt es jeweils eine kurze Pause. \n\n ' ...
    'Vor dem Block erscheint jeweils der Hinweis, ob Sie die Bilder in diesem Block' ...
    ' "BETRACHTEN" sollten ohne Ihre Gefuehle zu veraendern, oder die Geschichte der Bilder' ...
    ' "NEU ERZAEHLEN" sollten, damit sich das hervorgerufene Gefuehl abschwaecht. ' ...
    'Auf dem Bild erscheint zur Erinnerung ein "B" (betrachten) oder "N" (neu erzaehlen). ' ...
    'Im Anschluss bewerten Sie ihr aktuelles Empfinden mittels der folgenden Skala: ' ...
    '\n\n \n\n \n\n \n\n Mit dem linken Knopf koennen Sie auf der Skala nach links ' ...
    'navigieren, mit dem rechten Knopf nach rechts. Wenn Sie Ihre Bewertung abgeben wollen, ' ...
    'bestaetigen Sie bitte mit dem obersten Knopf auf der Button Box. \n\n \n\n ' ...
    'Wenn Sie soweit sind, druecken Sie den obersten Knopf auf der Box um die Aufgabe zu starten.'];


%% LOAD PICTURES and NAMES
folder = [pwd '/iaps_images']; % folder with all IAPS images
folder_scr = [pwd '/iaps_images_scr']; % folder with all scrambled IAPS images
folder_sam = [pwd '/sam_scale']; % folder with self-assessment manikin scales

% load the design (in which order images are presented to this participant,
% generated with esc1_er_allot.m)
er = load(['imagelist_er_' subject_id '.mat']); 
image_order = er.imgvect.blocks; 
trialtype = er.imgvect.type;
jitter = er.imgvect.jitter; % for timing the inter trial interval (ITI)


%% PREPARE DISPLAY
[windowHandle, rectangle] = Screen('OpenWindow', max(Screen('Screens')), ...
    [0.5 0.5 0.5]); 


% cell for loading the IAPS images
ratingImages = cell(nrRatings, 1);

for currRatingImage = 1 : nrRatings
    ratingImages{currRatingImage, 1} = double(imread([folder filesep ...
        num2str(image_order(currRatingImage)) '_mg.000_origDim.bmp' ]));
end

% cell for loading the scrambled images
ratingImagesScr = cell(nrRatings, 1);

for currRatingImage = 1 : nrRatings
    ratingImagesScr{currRatingImage, 1} = double(imread([folder_scr filesep ...
        num2str(image_order(currRatingImage)) '_mg.180_origDim.bmp'   ])); 
end

% cell for the SAM scale
ratingImagesSAM = cell(nrRatings,1);
if strcmp(scale_direction, 'l_to_r')
    for currRatingImage = 1 : nrRatings
        ratingImagesSAM{currRatingImage, 1} = double(imread([folder_sam filesep ...
            'SAM_L_to_R.bmp'   ]));
    end
end
if strcmp(scale_direction, 'r_to_l')
    for currRatingImage = 1 : nrRatings
        ratingImagesSAM{currRatingImage, 1} = double(imread([folder_sam filesep ...
            'SAM_R_to_L.bmp'   ]));
    end
end


%% PREPARE KEYBOARD
KbName('UnifyKeyNames');
escape_key = KbName('ESCAPE');

to_left = KbName('g'); % green button on diamond button box
to_right = KbName('b'); % blue button on diamond button box
% define a key to confirm the rating so that we can be sure the position
% was intentional, and not just the last cursor position before the rating
% period timed out:
choice_made = {'r'}; % red is the top button on the diamond button box

trigger_name = {'t'}; % signal the scanner sends for the trigger

choice_key = KbName(choice_made);
trigger_key = KbName(trigger_name);
escape_key = KbName('ESCAPE');

condition_cue = {'NEU ERZAEHLEN', 'BETRACHTEN'}; % reappraise or view condition
condition_cue_short = {'N', 'B'}; % short reminder for reappraise or view condition


%% PREPARE OUTPUT VARIABLES

timing = struct; %timing variable
responses = struct;
responses.responses = NaN(nrRatings,1); % responses for emotion ratings
% starting position for the rating: for sanity check - did participant move
% the cursor at all, or did they just confirm the starting position?
responses.starting_pos = NaN(nrRatings,1); 

timing.scramble_start_times = NaN(1,length(1:nrRatings)); % for pupil adaptation: record start time of the cue
timing.content_start_times = NaN(1,length(1:nrRatings));  % record start time for viewing or regulating emotion content
timing.trial_start_times = NaN(1,length(1:nrRatings)); % record start time for emotion ratings
timing.decision_times = NaN(1,length(1:nrRatings)); % time the rating was logged
timing.reaction_times = NaN(1,length(1:nrRatings)); % reaction time for the rating
% to make sure that all communication with the eyetracker 
% doesn't mess up the scanner timing:
timing.loopTimer = NaN(1,length(1:nrRatings)); 

timing.blank_start_times = NaN(nrRatings, 1); % record start times for ITI

timing.break_start_time = NaN(1,length(1:nrRatings)); % record start times for the breaks
timing.break_end_time = NaN(1,length(1:nrRatings)); % record end times for the breaks

timing.cue_start_time = NaN(1,length(1:nrRatings)); % record start times for the blocks (depiction of verbal cue before block)


%% PREPARE SCREEN SETTINGS

% determine center of the opened screen
[width, height] = Screen('WindowSize', max(Screen('Screens')));
centerX = floor(0.5 * width);
centerY = floor(0.5 * height);

% colors and text
txt_color = [200 200 200]; % off white
% fix_color = [127 127 127]; % middle gray
fix_color = [255 255 255]; % on top of the stimulus
box_color = [255 255 0]; % yellow
confirm_color = [0 255 255]; % turquoise

bg_color = GrayIndex(windowHandle, 0.5); % 0.5 = middle gray, 0 = black, 1 = white
wrapat_length = 60;
Screen('TextFont', windowHandle,'Arial');

txt_size_prompt = 18;
txt_size_task = 14;


% image height and width
% IAPS images:
imgH = 298;
imgW = 398;

% SAM images: 700 x 101 pixels
imgWS = 700;
imgHS = 101;

% selection frames on the SAM
frame1 = [rectangle(3)*1.91/8 rectangle(4)*2.7/6 rectangle(3)*2.36/8 rectangle(4)*3.3/6]; % frame position 1
frame2 = [rectangle(3)*2.38/8 rectangle(4)*2.7/6 rectangle(3)*2.82/8 rectangle(4)*3.3/6]; % frame position 2
frame3 = [rectangle(3)*2.84/8 rectangle(4)*2.7/6 rectangle(3)*3.29/8 rectangle(4)*3.3/6]; % frame position 3
frame4 = [rectangle(3)*3.31/8 rectangle(4)*2.7/6 rectangle(3)*3.76/8 rectangle(4)*3.3/6]; % frame position 4
frame5 = [rectangle(3)*3.78/8 rectangle(4)*2.7/6 rectangle(3)*4.22/8 rectangle(4)*3.3/6]; % frame position 5
frame6 = [rectangle(3)*4.24/8 rectangle(4)*2.7/6 rectangle(3)*4.69/8 rectangle(4)*3.3/6]; % frame position 6
frame7 = [rectangle(3)*4.71/8 rectangle(4)*2.7/6 rectangle(3)*5.15/8 rectangle(4)*3.3/6]; % frame position 7
frame8 = [rectangle(3)*5.17/8 rectangle(4)*2.7/6 rectangle(3)*5.62/8 rectangle(4)*3.3/6]; % frame position 8
frame9 = [rectangle(3)*5.64/8 rectangle(4)*2.7/6 rectangle(3)*6.08/8 rectangle(4)*3.3/6]; % frame position 9

frames = {frame1,frame2,frame3,frame4,frame5,frame6,frame7,frame8,frame9};


%% CREATE TEXTURES

% pics to display
ratingTexture = cell(nrRatings, 1);
for currRatingImage = 1 : nrRatings
    ratingTexture{currRatingImage, 1} = Screen('MakeTexture', ...
        windowHandle, ratingImages{currRatingImage, 1});
end

% pics to display scrambled
ratingTextureScr = cell(nrRatings, 1);
for currRatingImage = 1 : nrRatings
    ratingTextureScr{currRatingImage, 1} = Screen('MakeTexture', ...
        windowHandle, ratingImagesScr{currRatingImage, 1});
end

% sam scale to display
ratingTextureSAM = cell(nrRatings, 1);
for currRatingImage = 1 : nrRatings
    ratingTextureSAM{currRatingImage, 1} = Screen('MakeTexture', ...
        windowHandle, ratingImagesSAM{currRatingImage, 1});
end


%% DISPLAY PICTURES FOR RATING

try
    
    % drawing the instructions
    Screen('TextSize', windowHandle, txt_size_prompt);
    for i = 1:length(instructions) 
        
        Screen(windowHandle, 'FillRect', bg_color);
        DrawFormattedText(windowHandle, instructions{i}, 'center', 'center', txt_color, wrapat_length, [], [], 1.5);
        Screen('DrawTexture', windowHandle, ratingTextureSAM{currRatingImage, 1}, [], [centerX-350, centerY, centerX+350, centerY+101]);
        Screen(windowHandle, 'Flip');
        WaitSecs(.5);
        
        while 1
            [key_is_down, ~, key_code] = KbCheck;
            
            % flip the instructions with the "confirm" key!
            if key_is_down && any(key_code(choice_key)) 
                break;
            end
            
        end
        
    end
    
    Screen(windowHandle, 'FillRect', bg_color);
    
    DrawFormattedText(windowHandle, 'Warten auf Scanner', 'center', 'center', txt_color);
    Screen(windowHandle, 'Flip');
    
    FlushEvents;
    
    % wait for scanner trigger (expects a "t" input!)
    waiting_for_scanner = 1;
    
    while waiting_for_scanner
        
        [key_is_down, ~, key_code] = KbCheck;
        
        if key_is_down && any(key_code(trigger_key))
            
            waiting_for_scanner = 0;
            timing.start_time = GetSecs; % start of the scanner timing!
            
            if useEyetracker
                % send message to EyeLink recording file
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
        
        % first show cue for this block
        if currRatingImage == 1 || currRatingImage == 21 || currRatingImage == 41 || currRatingImage == 61 || currRatingImage == 81
            Screen(windowHandle, 'FillRect', bg_color);
            DrawFormattedText(windowHandle, condition_cue{trialtype(currRatingImage)}, 'center', 'center' , txt_color, wrapat_length, [], [], 1.5);
            Screen('flip',windowHandle);
            timing.cue_start_time(currRatingImage) = GetSecs;
            WaitSecs(1);
        end
        
        
        % DISPLAY SCRAMBLED IMAGE for adapting pupil response
        Screen('DrawTexture', windowHandle, ratingTextureScr{currRatingImage, 1});
        % with condition cue on top of the picture
        DrawFormattedText(windowHandle, condition_cue_short{trialtype(currRatingImage)}, centerX-6, centerY-11, txt_color, wrapat_length, [], [], 1.5);
        Screen('flip',windowHandle);
        
        loopTimer = tic();
        timing.scramble_start_times(currRatingImage) = GetSecs();
        WaitSecs(1); % 1 second to allow the pupil to adapt
        
        
        % REGULATION / VIEW SCREEN
        Screen('TextSize', windowHandle, txt_size_task);
        Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1});
        % with condition cue on top of the picture
        DrawFormattedText(windowHandle, condition_cue_short{trialtype(currRatingImage)}, centerX-6, centerY-11, txt_color, wrapat_length, [], [], 1.5);
        
        Screen('flip',windowHandle);
        timing.content_start_times(currRatingImage) = GetSecs;
        WaitSecs(7); % 7 seconds content display
        
        
        % EMOTION RATINGS   
        Screen('TextSize', windowHandle, txt_size_task);
        Screen('DrawTexture', windowHandle, ratingTextureSAM{currRatingImage, 1});
        % with condition cue above the picture
        DrawFormattedText(windowHandle, 'AKTUELLES GEFUEHL', 'center', centerY-70 , txt_color, wrapat_length, [], [], 1.5);
        
        % on top of the scale: draw yellow square that randomly pops up in
        % one spot
        start_rating = randi(9);
        responses.starting_pos(currRatingImage) = start_rating;
        Screen('FrameRect', windowHandle, box_color, frames{start_rating}, 5);
        Screen('flip',windowHandle);
        timing.trial_start_times(currRatingImage) = GetSecs;
        
        % move the selection frame left or right to answer, confirm the
        % answer with red button 'r'
        
        % needs to be outside while loop! otherwise cursor does not move
        pointerpos_curr = start_rating; 
        
        % read out answer
        choice_made = false;
        flipscreen = false;
        while ~choice_made
            
            [key_is_down, ~, key_code] = KbCheck;
            
            if any(key_code(to_left)) && pointerpos_curr > 1; % && ensures that cursor stops at 0
                
                % arrow left: shift cursor with left arrow to next frame to the
                % left until space was pressed or left end of the scale was reached
                Screen('DrawTexture', windowHandle, ratingTextureSAM{currRatingImage, 1});
                DrawFormattedText(windowHandle, 'AKTUELLES GEFUEHL', 'center', centerY-70, txt_color, wrapat_length, [], [], 1.5);
                pointerpos_curr = pointerpos_curr - 1;
                pointer_move_l = pointerpos_curr;
                Screen('FrameRect', windowHandle, box_color, frames{pointer_move_l}, 5);
                Screen('flip',windowHandle);
                WaitSecs(0.2);
                
            elseif any(key_code(to_right)) && pointerpos_curr < 9;
                % shift cursor with right arrow to next separator 
                % to the right until space is pressed or position
                % 9 was reached
                Screen('DrawTexture', windowHandle, ratingTextureSAM{currRatingImage, 1});
                DrawFormattedText(windowHandle, 'AKTUELLES GEFUEHL', 'center', centerY-70 , txt_color, wrapat_length, [], [], 1.5);
                pointerpos_curr = pointerpos_curr + 1;
                pointer_move_r = pointerpos_curr;
                Screen('FrameRect', windowHandle, box_color, frames{pointer_move_r}, 5);
                Screen('flip',windowHandle);
                WaitSecs(0.2);
                
            elseif key_is_down && any(key_code(choice_key))&& length(find(key_code))<2;
                % confirm answer by pressing the red button on the button box
                if strcmp(scale_direction, 'l_to_r') 
                    % sad left, happy right on the 9 point SAM scale
                    responses.responses(currRatingImage) = pointerpos_curr;
                end
                if strcmp(scale_direction, 'r_to_l') 
                    % happy left, sad right on the 9 point SAM scale
                    responses.responses(currRatingImage) = 10-pointerpos_curr;
                end
                timing.decision_times(currRatingImage) = GetSecs;
                timing.reaction_times(currRatingImage) = timing.decision_times(currRatingImage) - timing.trial_start_times(currRatingImage);
                choice_made = true;
                break;
            end
            
            if key_is_down && key_code(escape_key)
                save(results_file_name, 'responses', 'timing', 'jitter');
                break;
            end
            
            if (GetSecs - timing.trial_start_times(currRatingImage)) >= rating_time % 4 seconds rating time
                disp('no rating made')
                break;
            end
        end %while loop choice
        
        if ~isnan(responses.responses(currRatingImage))
            % display answer feedback if an answer was made:
            % selection frame changes its color to yellow
            Screen('DrawTexture', windowHandle, ratingTextureSAM{currRatingImage, 1});
            DrawFormattedText(windowHandle, 'AKTUELLES GEFUEHL', 'center', centerY-70 , txt_color, wrapat_length, [], [], 1.5);
            Screen('FrameRect', windowHandle, confirm_color, frames{pointerpos_curr}, 5);
            Screen('flip',windowHandle);
            WaitSecs(0.1); % display time answer feedback % this "wait" must
            % live inside this IF statement, otherwise getting too slow 
            % if the participant does not answer!
        end
        
        %fixation cross - jittered ITI
        Screen('TextSize', windowHandle, txt_size_task);
        Screen(windowHandle, 'FillRect', bg_color);
        %DrawFormattedText(windowHandle, '+', 'center', 'center', txt_color);
        Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
        Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
        
        Screen('flip', windowHandle);
        timing.blank_start_times(currRatingImage) = GetSecs;
        
        
        save(results_file_name, 'responses',  'timing'); %save here to update results on every loop
   
        
        if useEyetracker
            % send message to EyeLink recording file
            Eyelink('message', sprintf('trial: %d', currRatingImage));
            Eyelink('command', sprintf('record_status_message ''End trial %d''', currRatingImage));
        end
        
        % now time the whole loop - it should not be longer, including the
        % jitter, than we planned it to be for the total trial:
        % cue present time + pic present time + rating time + jitter
        while toc(loopTimer) < 1 + 7 + 4 + jitter(currRatingImage), end
        timing.loopTimer(currRatingImage) = toc(loopTimer);
        
        
        % after each block, allow for small break of 15 seconds
        if currRatingImage == 20 || currRatingImage == 40 || currRatingImage == 60 || currRatingImage == 80
            timer2 = tic();
            timing.break_start_time(currRatingImage) = GetSecs;
            
            if useEyetracker == true
                %send message to EyeLink recording file
                Eyelink('message', 'StartPause');
                %Eyelink('message', 'CUE_%d', event);
                Eyelink('command', 'record_status_message ''StartPause''');
            end
            
            Screen('TextSize', windowHandle, txt_size_task);
            Screen(windowHandle, 'FillRect', bg_color);
            counter = break_time;
            
            %%
            while toc(timer2) <= break_time
                DrawFormattedText(windowHandle, 'PAUSE', 'center', centerY-70, txt_color);
                DrawFormattedText(windowHandle, num2str(break_time - ceil(toc(timer2))), 'center', 'center', txt_color);
                
                Screen('flip', windowHandle);
                % %                 WaitSecs(1);
                counter = counter-1;
            end
            
            %%
            timing.break_end_time(currRatingImage) = GetSecs;
            
        end
        
    end % end of main draw loop
    
    if useEyetracker == true
        %send message to EyeLink recording file
        Eyelink('message', 'ExpEnd');
        Eyelink('command', 'record_status_message ''ExpEnd''');
    end
    
    %fixation at the end of the run
    Screen('TextSize', windowHandle, txt_size_task);
    Screen(windowHandle, 'FillRect', bg_color);
    Screen('DrawLine', windowHandle, txt_color, centerX-5, centerY, centerX+5, centerY, 2);
    Screen('DrawLine', windowHandle, txt_color, centerX, centerY-5, centerX, centerY+5, 2);
    Screen('flip', windowHandle);
    WaitSecs(10*TR);

    save(results_file_name, 'responses',  'timing'); 
   
catch
    
    psychrethrow(psychlasterror);
    
end %end of try-catch


%% CLOSE DISPLAY
Screen('CloseAll');

end
