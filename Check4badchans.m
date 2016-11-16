clear all
close all
clc

path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\PortableTestJames\'; %Monash PC

subject_folder = {'James'};
allsubj = {'PortableJamesTest'};

Other_System = {}; %List participants tested using a diferent EEG system

Monash_System = {'James'};%List participants tested using the 32 channel 
                           %BrainProducts’ portable LiveAmp EEG system

duds = [];
single_participants = [1];

file_start = 1;

blocks = {[1:9]};


if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    blocks([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    blocks = blocks(single_participants);
end

nchan = 35;

% for s=file_start:length(subject_folder)
%     matfiles{s} = [path_temp subject_folder{s} '\' allsubj{s} 'chanvars.mat'];
% end
for s=file_start:length(subject_folder)
    matfiles{s} = [path_temp subject_folder{s} '\' allsubj{s} 'chanvars.mat'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=file_start:length(subject_folder)
    if ismember(subject_folder{s},Other_System)
        chanlocs = readlocs('cap64.loc');
    else
       chanlocs = readlocs ('actiCAP32_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
    end
    disp(subject_folder{s})
    load(matfiles{s})

    chanVar = double(chanVar);
    
    badchans = [];
    changechans = []; % must be in same order as badchans.
    chanVar(badchans(1:end),:) = chanVar(changechans(1:end),:);
    
    avVar = mean(chanVar,2); 
        
    % average variance for each channel across all 16 conditions
    % on a second sweep for a given subject, might want to plot topo again
    % after getting rid of a really bad one (to make it easier to see other
    % bad channels) - so do something like:
    % avVar(104) = avVar(103);  % quick hack - make a reall bad chan equal its neighbor
    
    figure;
    topoplot(avVar,chanlocs,'plotchans',[1:32],'electrodes','numbers');
    title(subject_folder{s})
    
    figure; hold on
    h = plot(chanVar(1:32,:));
    title(subject_folder{s})
    legend(h,'Location','NorthEast');
    
end