function [responses] = ratings_escf1(scale_direction, results_file_name, ...
    rating_type, language)
% health and taste ratings for the food images in the ESC fMRI study

% inputs: 
% scale_direction: 'l_to_r' for left to right or 'r_to_l' for right to left
% (negative -> positive directionality)
% results_file_name: allocated in master function that called this script
% (results_health_[subID].mat or results_taste_[subID].mat)
% rating_type: as a character string: 'Gesundheit' (health) or 'Geschmack'
% (taste) - careful with the spelling, must match language switch below!
% language: as a character string: 'DE' (for German), 'EN' (for English)


%% PREPARE INSTRUCTIONS

switch language
    
    case 'DE'
        if strcmp(rating_type, 'Appetit') == 1
            instructions{1} = 'BEWERTEN: \n\n \n\n Ab jetzt gehen Ihre Bewertungen in die Studie ein. In dieser Runde bitten wir Sie anzugeben wieviel APPETIT Sie auf die gezeigten Nahrungsmittel haben. Das heisst, Sie beantworten fuer sich die Frage: "Wie gern wuerde ich dieses Nahrungsmittel am Ende der Studie essen?". \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit den Instruktionen fortzufahren.';
            
            instructions{2} = 'Sie werden auf dem Bildschirm jeweils 1 Nahrungsmittel sehen und darunter eine Bewertungsskala. \n\n Wichtig: Sie muessen eine Maustaste gedrueckt halten, damit die rote Linie auf der Skala erscheint, die Sie mit der Maus hin und her bewegen koennen. Diese Linie zeigt Ihre Bewertung an. Fuehren Sie die rote Linie zu der Stelle auf der Bewertungsskala, die ihrer Bewertung am besten entspricht. Bitte bestaetigen Sie ihre Bewertung dann mit der Leertaste. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit den Instruktionen fortzufahren.';
            
            instructions{3} = 'Sollten Sie sich nicht ganz sicher sein und Ihre Bewertung nicht innerhalb von 3 Sekunden abgeben wird das Bild am Ende der Runde noch einmal wiederholt. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit den Instruktionen fortzufahren.';
            
            instructions{4} = 'Auf dem naechsten Bildschirm startet die Bewertungsaufgabe. Ueber der Skala sehen Sie immer angegeben, was Sie gerade bewerten. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit der Aufgabe zu beginnen.';
        end
        
        if strcmp(rating_type, 'Gesundheit') == 1
            instructions{1} = 'BEWERTEN: \n\n \n\n In diesem Teil werden Sie beurteilen, wie GESUND Sie verschiedene Nahrungsmittel finden. Der Geschmack soll dabei keine Rolle spielen. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit den Instruktionen fortzufahren.';
            
            instructions{2} = 'Bitte beurteilen Sie die Gesundheit der Nahrungsmittel nach Ihrem eigenen Empfinden. Es gibt keine richtigen oder falschen Antworten. Als Portion gilt die Menge, die Ihrer geballten Faust entspricht, bei Schokolade und Keksen ist es eine halbe Faust. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit der Aufgabe zu beginnen.';
        end
        
        if strcmp(rating_type, 'Geschmack') == 1
            instructions{1} = 'BEWERTEN: \n\n \n\n In diesem Teil werden Sie beurteilen, wie gut Ihnen verschiedene Nahrungsmittel SCHMECKEN. Dabei soll es keine Rolle spielen, wie gesund Sie das Nahrungsmittel finden. Bitte geben Sie kein generelles Urteil ab, sondern beurteilen Sie das Produkt, wie Sie es VOR SICH SEHEN und wie es Ihnen JETZT schmecken wuerde, wenn Sie es im Anschluss an den Versuch zu essen bekaemen. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit den Instruktionen fortzufahren.';
            
            instructions{2} = 'Falls Sie ein bestimmtes Produkt nicht kennen geben Sie bitte Ihre beste Schaetzung an, ob Sie das Nahrungsmittel moegen wuerden. \n\n Bitte antworten Sie nach Ihrem Empfinden. Es gibt keine richtigen oder falschen Antworten. \n\n \n\n Bitte druecken Sie die rechte Pfeiltaste um mit der Aufgabe zu beginnen.';
        end
        
    case 'EN'
        if strcmp(rating_type, 'Health') == 1
            instructions{1} = 'RATINGS: \n\n \n\n In this part you will rate how HEALTHY certain food items are, regardless of their taste. \n\n \n\n Please press the right arrow on the keyboard to continue.';
            
            instructions{2} = 'Please give us your honest opinion. There are no right or wrong answers. The portion size is one handful of the item, for chocolate and cookies it is half a handful. \n\n \n\n Please press the right arrow on the keyboard to start the task.';
        end
        
        if strcmp(rating_type, 'Taste') == 1
            instructions{1} = 'RATINGS: \n\n \n\n In this part you will rate how much you like the TASTE of certain food items, regardless of their healthiness. Please do not give us a general impression. Instead, please indicate how tasty you would find the item that you see ON THE SCREEN if you got to eat it TODAY after the session. \n\n \n\n Please press the right arrow on the keyboard to continue.';
            
            instructions{2} = 'If you are unsure about a specific product or brand, please give us your best guess how much you would like the taste of this item. \n\n Please give us your honest opinion. There are no right or wrong answers. \n\n \n\n Please press the right arrow on the keyboard to start the task.';
        end
        
end


%% LOAD PICTURES

folder = [pwd '/food_images'];
filenames = dir([folder '/*.bmp']);
nrRatings = size(filenames, 1);
% shuffle images so that everyone sees them in a different order
image_order = randperm(nrRatings); 
ratingImages = cell(nrRatings, 1);

for currRatingImage = 1 : nrRatings
    ratingImages{currRatingImage, 1} = double(imread([folder '/' ...
        filenames(image_order(currRatingImage)).name]));
end


%% PREPARE DISPLAY

[windowHandle, rectangle] = Screen('OpenWindow', max(Screen('Screens')), [0 0 0], [0 0 1280 768]);

% prepare cell with textures for the food images to rate
ratingTexture = cell(nrRatings, 1);
for currRatingImage = 1 : nrRatings
    ratingTexture{currRatingImage, 1} = Screen('MakeTexture', ...
        windowHandle, ratingImages{currRatingImage, 1});
end


%% PREPARE KEYBOARD

KbName('UnifyKeyNames');
escape_key = KbName('ESCAPE');
choice = KbName('space');
moveon = KbName('RightArrow');


%% PREPARE OUTPUT VARIABLES

timing = struct;
responses = NaN(nrRatings,1);


%% PREPARE SCREEN SETTINGS

txt_color = [255 255 255]; %white
pointer_color = [255 0 0]; %red
bg_color = BlackIndex(windowHandle);
% txt_size_prompt = 20;
wrapat_length = 60;
Screen('TextSize', windowHandle, 20);
Screen('TextFont', windowHandle,'Arial');


%% DISPLAY PICTURES FOR RATING

timing.start_time = GetSecs;
timing.trial_start_times = 1:nrRatings;
timing.decision_times = 1:nrRatings;
timing.blank_start_times = 1:nrRatings;
timing.reaction_times = 1:nrRatings;

for i=1:length(instructions) %drawing the instructions
    Screen(windowHandle, 'FillRect', bg_color);
    DrawFormattedText(windowHandle, instructions{i}, 'center', 'center', txt_color, wrapat_length, [], [], 1.5);
    Screen(windowHandle, 'Flip');
    WaitSecs(.5);
    while 1
        [key_is_down, ~, key_code] = KbCheck;
        if key_is_down && any(key_code(moveon))
            break;
        end
    end
end
Screen(windowHandle, 'FillRect', bg_color);
FlushEvents;


for currRatingImage = 1 : nrRatings % main draw loop
    
    % currRatingImage = image_order(currRatingImage); % to re-sort ratings in the order of original image numbers (1:180 = a-z)
    
    % image
    Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1});
    
    % slider bar
    drawscale_escf1(rectangle, windowHandle, txt_color, scale_direction)
    % explanation which rating is collected
    barPos = [rectangle(3)*2/6 rectangle(4)*1.02/1.11 rectangle(3)*4/6 rectangle(4)*1.02/1.11];
    txtPos = barPos(2)-100;
    DrawFormattedText(windowHandle, rating_type, 'center', txtPos, txt_color, wrapat_length, [], [], 1.5);
    
    % mouse shows up somewhere random on the screen
    SetMouse((randInt(1, 1, (round(barPos(1)):round(barPos(3))), 'noreplace')), (randInt(1, 1, (round(rectangle(4)*1.02/1.25):round(rectangle(4)*1.02/1.18)), 'noreplace')), max(Screen('Screens')));
    
    Screen('flip', windowHandle);
    
    % timing
    timing.trial_start_times(currRatingImage) = GetSecs;
    
    
    % update the position of the red pointer on the slider bar until a
    % choice is confirmed by pressing space
    choice_made = false;
    
    button_was_pressed = false;
    
    while ~choice_made
        
        % read out the mouse only within the boundaries of the slider bar
        [mouseXpos, ~, buttons] = GetMouse(max(Screen('Screens')));
        if mouseXpos < round(rectangle(3)*2/6) % left bound
            mouseXpos = round(rectangle(3)*2/6);
        end
        if mouseXpos > round(rectangle(3)*4/6)% right bound
            mouseXpos = round(rectangle(3)*4/6);
        end
        
        Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1});
        
        
        if any(buttons) % scale is only drawn if mouse button is held down
            button_was_pressed = true;
            drawscale_escf1(rectangle, windowHandle, txt_color, scale_direction)
            barPos = [rectangle(3)*2/6 rectangle(4)*1.02/1.11 rectangle(3)*4/6 rectangle(4)*1.02/1.11];
            txtPos = barPos(2)-100;
            DrawFormattedText(windowHandle, rating_type, 'center', txtPos, txt_color, wrapat_length, [], [], 1.5);
            Screen('DrawLine', windowHandle, pointer_color, mouseXpos, rectangle(4)*1.04/1.11, mouseXpos, rectangle(4)/1.15, 3);
            Screen('flip',windowHandle);
        end
                
        [key_is_down, ~, key_code, ~] = KbCheck;
        
        if (key_is_down && key_code(choice)) && button_was_pressed % && any(buttons) %space is pressed to confirm choice, sanity check whether subject did move the mouse at all
            switch scale_direction
                case 'l_to_r'
                    responses(currRatingImage) = round(mouseXpos) - round(rectangle(3)*2/6); % current x - left x
                case 'r_to_l'
                    responses(currRatingImage) = abs( (round(mouseXpos) - round(rectangle(3)*2/6)) - (round(rectangle(3)*4/6) - round(rectangle(3)*2/6)) ); % subtract length of rating bar
            end
            timing.decision_times(currRatingImage) = GetSecs;
            timing.reaction_times(currRatingImage) = timing.decision_times(currRatingImage) - timing.trial_start_times(currRatingImage);
            choice_made = true;
            break;
        end
        
        
        if (GetSecs - timing.trial_start_times(currRatingImage)) >=3 %3 seconds picture presentation time
            break;
        end
        
    end % end of while loop
    
    %fixation cross
    Screen(windowHandle, 'FillRect', bg_color);
    Screen('DrawLine', windowHandle ,txt_color, rectangle(3)/2, rectangle(4)*4.25/9, rectangle(3)/2, rectangle(4)*4.75/9, 3);
    Screen('DrawLine', windowHandle ,txt_color, rectangle(4)*4.25/9+(rectangle(3)-rectangle(4))/2, rectangle(4)/2, rectangle(4)*4.75/9+(rectangle(3)-rectangle(4))/2, rectangle(4)/2, 3);
    
    Screen('flip', windowHandle);
    timing.blank_start_times(currRatingImage) = GetSecs;
    
    WaitSecs(1.5);  % time spent in between ratings
    
    save(results_file_name, 'responses', 'timing', 'image_order', 'scale_direction'); %save here to update results on every loop
    
    
end % end of main draw loop

% check whether pictures were missed
no_responses = isnan(responses);
timing.decision_reps=ones(nrRatings,1);
timing.decision_reps=timing.decision_reps+no_responses;
missing=1:nrRatings;


% draw loop for missed pictures
while sum(no_responses) > 0
    
    for currRatingImage = missing(no_responses)
                
        % image
        Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1});
        
        % slider bar
        drawscale_escf1(rectangle, windowHandle, txt_color, scale_direction)
        % explanation which rating is collected
        barPos = [rectangle(3)*2/6 rectangle(4)*1.02/1.11 rectangle(3)*4/6 rectangle(4)*1.02/1.11];
        txtPos = barPos(2)-100;
        DrawFormattedText(windowHandle, rating_type, 'center', txtPos, txt_color, wrapat_length, [], [], 1.5);
        
        % mouse shows up somewhere random on the screen
        SetMouse((randInt(1, 1, (round(barPos(1)):round(barPos(3))), 'noreplace')), (randInt(1, 1, (round(rectangle(4)*1.02/1.25):round(rectangle(4)*1.02/1.18)), 'noreplace')), max(Screen('Screens')));
        
        Screen('flip', windowHandle);
        
        % timing
        timing.trial_start_times(currRatingImage) = GetSecs;
        
        
        % update the position of the red pointer on the slider bar until a
        % choice is confirmed by pressing space
        choice_made = false;
        
        button_was_pressed = false;
        
        while ~choice_made
            
            % read out the mouse only within the boundaries of the slider bar
            [mouseXpos, ~, buttons] = GetMouse(max(Screen('Screens')));
            if mouseXpos < round(rectangle(3)*2/6) % left bound
                mouseXpos = round(rectangle(3)*2/6);
            end
            if mouseXpos > round(rectangle(3)*4/6)% right bound
                mouseXpos = round(rectangle(3)*4/6);
            end
            
            Screen('DrawTexture', windowHandle, ratingTexture{currRatingImage, 1});
            
            
            if any(buttons) % scale is only drawn if mouse button is held down
                button_was_pressed = true;
                drawscale_escf1(rectangle, windowHandle, txt_color, scale_direction)
                barPos = [rectangle(3)*2/6 rectangle(4)*1.02/1.11 rectangle(3)*4/6 rectangle(4)*1.02/1.11];
                txtPos = barPos(2)-100;
                DrawFormattedText(windowHandle, rating_type, 'center', txtPos, txt_color, wrapat_length, [], [], 1.5);
                Screen('DrawLine', windowHandle, pointer_color, mouseXpos, rectangle(4)*1.04/1.11, mouseXpos, rectangle(4)/1.15, 3);
                Screen('flip',windowHandle);
            end
                        
            [key_is_down, ~, key_code, ~] = KbCheck;
            
            if (key_is_down && key_code(choice)) && button_was_pressed % && any(buttons) % space is pressed to confirm choice
                switch scale_direction
                    case 'l_to_r'
                        responses(currRatingImage) = round(mouseXpos) - round(rectangle(3)*2/6); % current x - left x
                    case 'r_to_l'
                        responses(currRatingImage) = abs( (round(mouseXpos) - round(rectangle(3)*2/6)) - (round(rectangle(3)*4/6) - round(rectangle(3)*2/6)) ); % subtract length of rating bar
                end
                timing.decision_times(currRatingImage) = GetSecs;
                timing.reaction_times(currRatingImage) = timing.decision_times(currRatingImage) - timing.trial_start_times(currRatingImage);
                choice_made = true;
                break;
            end
                     
            if (GetSecs - timing.trial_start_times(currRatingImage)) >=4 %3 seconds picture presentation time
                break;
            end
            
        end % end of while loop
        
        %fixation cross
        Screen(windowHandle, 'FillRect', bg_color);
        Screen('DrawLine', windowHandle ,txt_color, rectangle(3)/2, rectangle(4)*4.25/9, rectangle(3)/2, rectangle(4)*4.75/9, 3);
        Screen('DrawLine', windowHandle ,txt_color, rectangle(4)*4.25/9+(rectangle(3)-rectangle(4))/2, rectangle(4)/2, rectangle(4)*4.75/9+(rectangle(3)-rectangle(4))/2, rectangle(4)/2, 3);

        Screen('flip', windowHandle);        
        timing.blank_start_times(currRatingImage) = GetSecs;
        
        save(results_file_name, 'responses', 'timing', 'image_order', 'scale_direction'); %save here to update results on every loop
        
        WaitSecs(1.5);  % time spent in between pictures
        
        
    end % end no response draw loop
    
    no_responses = isnan(responses);
    xtemp= isnan(responses);
    timing.decision_reps=timing.decision_reps+xtemp; % counter for number of times presented
end % catch no response while loop

% sort the ratings alphabetically (A:Z corresponds to image 1:180)
responses_in_image_order = responses;
responses = [];
responses(image_order) = responses_in_image_order;


save(results_file_name, 'responses', 'responses_in_image_order', 'timing', 'image_order', 'scale_direction');

%% CLOSE DISPLAY
Screen('CloseAll');

end
