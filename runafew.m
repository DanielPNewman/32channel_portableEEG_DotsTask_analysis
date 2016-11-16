%%% runafew
clear
close all
clc

%%
path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\PortableTestJames\'; %Monash PC

%%
subject_folder = {'James'}; 
allsubj = {'PortableJamesTest'};

%%

Other_System = {}; %List participants tested using a diferent EEG system
Monash_System = {'James'};%List participants tested using the 32 channel 
                           %BrainProducts’ portable LiveAmp EEG system
%%
duds = []; 
single_participants = [];

allblocks = {[1:9]};

allbadchans = {[21]}; % James. electrode 21 was turned off

%% Used the Eyelink eyetracker?
use_el=0; %1=yes; 0=no

%%
if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    allblocks([duds]) = [];
    allbadchans([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder([single_participants]);
    allsubj = allsubj([single_participants]);
    allblocks = allblocks([single_participants]);
    allbadchans = allbadchans([single_participants]);
end

%% CSD
if ~isempty(Other_System)
E = textread('chans64.asc','%s');
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
% MapMontage(M);
[G_TCD,H_TCD] = GetGH(M);
end

E = textread('chans32_monash.asc','%s');
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
% MapMontage(M);
[G_monash,H_monash] = GetGH(M);

for s=1:length(allsubj)
    disp(allsubj{s})
    
    blocks = allblocks{s};
    badchans = allbadchans{s};
    clear paths files matfiles ET_files ET_matfiles; k=0;
    for n=1:length(blocks)
        k=k+1;
        if ismember(subject_folder{s},Other_System)
            files{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '.bdf'];
            matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '.mat'];
            ET_files{k}=[path_temp 'Samples_and_Events\' allsubj{s} '_' num2str(blocks(n)) '.asc'];
            ET_matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '_ET.mat'];
        elseif ismember(subject_folder{s},Monash_System)
            files{k} = [allsubj{s} '_',num2str(blocks(n)) '.vhdr'];
            paths{k} = [path_temp subject_folder{s} '\'];
            matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '.mat'];
            ET_files{k}=[path_temp 'Samples_and_Events\' allsubj{s} '_' num2str(blocks(n)) '.asc'];
            ET_matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '_ET.mat'];
        end
    end            

    if ismember(subject_folder{s},Other_System)
        G_CSD = G_TCD;
        H_CSD = H_TCD;
        Other_System_preprocess
    elseif ismember(subject_folder{s},Monash_System)
        G_CSD = G_monash;
        H_CSD = H_monash;
        Monash_System_preprocess
    end
end
