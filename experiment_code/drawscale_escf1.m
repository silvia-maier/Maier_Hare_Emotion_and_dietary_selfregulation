function drawscale_escf1(rectangle, windowHandle, txt_color, scale_direction, language)
% draws the rating scale for the food taste and health ratings

% Inputs: 
% rectangle = screen settings (pixel width, height) from Psychtoolbox call
% that opened the window
% windowHandle = pointer to the screen (number)
% txt_color = color of text depiction (e.g. [0 0 0] = black)
% scale_direction = as a character string: 'l_to_r' for left to right scale
% (negative -> positive) or 'r_to_l' (for right to left order of negative -> positive)

% paint the scale in this rectangle
scaleframe = [rectangle(3)*2/6, rectangle(4)/1.13, rectangle(3)*4/6, rectangle(4)/1.11];

% another rectangle on top of the rating bar that goes from -0.5 to 0.5 to
% depict the neutral zone
neutralframe = [rectangle(3)*31.9/66, rectangle(4)/1.13, rectangle(3)*34.1/66, rectangle(4)/1.11]; 

Screen('FillRect', windowHandle, [127 127 127], scaleframe, 4)
Screen('FillRect', windowHandle, [200 200 0], neutralframe, 4)

% draw explanation of scale (zero in the middle)

switch scale_direction
    case 'l_to_r'
        [~, ~] = DrawFormattedText(windowHandle,'-5', (rectangle(3)*2/6)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-4', (rectangle(3)*24.2/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-3', (rectangle(3)*26.4/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-2', (rectangle(3)*28.6/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-1', (rectangle(3)*30.8/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'Neutral', (rectangle(3)*3/6)-30, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+1', (rectangle(3)*35.2/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+2', (rectangle(3)*37.4/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+3', (rectangle(3)*39.6/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+4', (rectangle(3)*41.8/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+5', (rectangle(3)*4/6)-10, rectangle(4)*1.02/1.1, txt_color);
    case 'r_to_l'
        [~, ~] = DrawFormattedText(windowHandle,'+5', (rectangle(3)*2/6)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+4', (rectangle(3)*24.2/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+3', (rectangle(3)*26.4/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+2', (rectangle(3)*28.6/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'+1', (rectangle(3)*30.8/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'Neutral', (rectangle(3)*3/6)-30, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-1', (rectangle(3)*35.2/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-2', (rectangle(3)*37.4/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-3', (rectangle(3)*39.6/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-4', (rectangle(3)*41.8/66)-10, rectangle(4)*1.02/1.1, txt_color);
        [~, ~] = DrawFormattedText(windowHandle,'-5', (rectangle(3)*4/6)-10, rectangle(4)*1.02/1.1, txt_color);
end


end