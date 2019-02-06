function esc1_er_allot(subject_id, rand_type)
% builds the design for the emotion paradigm

% Inputs:
% subject_id = subject identifier as a character string, e.g. '9999'
% rand_type = double (from 1 to 5) that indicates in which pseudorandom
% order (counterbalanced across men and women) the task should be performed
% this pseudorandomization also made sure that half the participants of both genders
% performed the emotion regulation task first and food choice second, or
% vice versa

% Outputs:
% imgvect.blocks with a 1:100 column vector for the IAPS images
% that should be presented
% imgvect.jitter with a 1:100 column vector for the jitter (1-5 seconds uniform)
% imgvect.type with a 1:100 column vector for the type of viewing condition
% - viewing conditions: 1 = regulate, 2 = view (always downregulate)

% define the output filename
filename = ['imagelist_er_' subject_id '.mat'];

% IAPS allocation / randomization

% identifier in the IAPS set:

% negative set A
negsetA = Shuffle([1300
    2055.1
    2095
    2981
    3015
    3181
    3301
    3550
    6020
    6212
    6370
    6540
    6838
    9040
    9180
    9181
    9265
    9435
    9520
    9570]);

% negative set B
negsetB = Shuffle([1525
    2352.2
    2683
    2800
    3051
    3250
    3530
    6312
    6415
    6570.1
    9140
    9250
    9252
    9253
    9300
    9430
    9561
    9571
    9635.1
    9800]);

% positive set A
possetA = Shuffle([1460
    1710
    1721
    1750
    1810
    1920
    2050
    2080
    2091
    2224
    2260
    2311
    2351
    2375.2
    2550
    5600
    5626
    5831
    5890
    8190]);

% positive set B
possetB = Shuffle([1440
    1463
    1720
    1731
    2058
    2071
    2150
    2170
    2303
    2345
    2620
    2655
    2660
    5390
    5594
    5628
    5830
    7580
    8461
    8497]);

% neutral set
neutrset = Shuffle([2020
    2200
    2214
    2357
    2480
    2493
    2570
    2880
    2890
    7030
    7036
    7150
    7161
    7170
    7186
    7224
    7590
    7705
    7830
    8030]);

% randomly allocate whether set A will be regulate or view (B vice versa)
randind = rand > 0.5;
if randind == 1
    posset1 = possetA; %regulate
    posset2 = possetB; %view
    
    negset1 = negsetA; %regulate
    negset2 = negsetB; %view
else
    posset1 = possetB; %regulate
    posset2 = possetA; %view
    
    negset1 = negsetB; %regulate
    negset2 = negsetA; %view
end

% allocate blocks
% trials: 1 sec cue + 7 sec regulate + 4 sec rating + 1-5 sec uniform jitter
% to make sure that block length is the same
% => 15 sec on average per trial

% put all trials into one block to obtain a 5 minute recording for emotion
% regulation type
% make sure that not two regulate blocks are following each other

% set 1 = regulate
% set 2 = view
% neutral is always view

switch rand_type 
    case 1 %neg view, pos reg, neutral, neg reg, pos view
        imgvect.blocks = cat(1, negset2, posset1, neutrset, negset1, posset2);
        type_tmp = repmat([2 1 2 1 2], 20, 1);
        imgvect.type = type_tmp(:);
    case 2 %pos view, neg reg, neutral, pos reg, neg view
        imgvect.blocks = cat(1, posset2, negset1, neutrset, posset1, negset2);
        type_tmp = repmat([2 1 2 1 2], 20, 1);
        imgvect.type = type_tmp(:);
    case 3 %pos reg, neg view, neutral, pos view, neg reg
        imgvect.blocks = cat(1, posset1, negset2, neutrset, posset2, negset1);
        type_tmp = repmat([1 2 2 2 1], 20, 1);
        imgvect.type = type_tmp(:);
    case 4 %neg reg, pos view, neutral, neg view, pos reg
        imgvect.blocks = cat(1, negset1, posset2, neutrset, negset2, posset1);
        type_tmp = repmat([1 2 2 2 1], 20, 1);
        imgvect.type = type_tmp(:);
    case 5 %pos view, neg reg, neutral, neg view, pos reg
        imgvect.blocks = cat(1, posset2, negset1, neutrset, negset2, posset1);
        type_tmp = repmat([2 1 2 2 1], 20, 1);
        imgvect.type = type_tmp(:);
end

% determine the jitter between trials (uniform distribution)
jitter_tmp = repmat(1:5,20,1);
imgvect.jitter = Shuffle(jitter_tmp(:));

save(filename, 'imgvect', 'rand_type', 'randind')

end




