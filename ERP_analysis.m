clear
close all
clc
chanlocs = readlocs('actiCAP32_ThetaPhi.elp','filetype','besa'); %DN for actiCAP

path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\PortableTestJames\'; %Monash PC

%% Create paths for subjects

subject_folder = {'James'};
allsubj = {'PortableJamesTest'};

%% Define if subjects collected using Monash system (the 32 channel BrainProducts’ portable LiveAmp EEG system)

Monash_System = {'James'};%List participants tested using the 32 channel
%BrainProducts’ portable LiveAmp EEG system

Other_System = {}; %List participants tested using a diferent EEG system

%% Was saliva collected and DAT1 genotypes extracted?
DAT1_genotypes=0; %1=yes; 0=no;

%% Load DAT1 genotypes if you collected them
if DAT1_genotypes
    load('DAT1genotypesforMatlab.mat')
    DAT1_split=[]; DAT1_nosplit=[]; dud_temp=[];
    for s = 1:length(subject_folder)
        for i = 1:size(DAT1genotypesforMatlab,1)
            if strcmp(subject_folder{s},DAT1genotypesforMatlab{i,1})
                if ~isempty(DAT1genotypesforMatlab{i,2})
                    DAT1_split(s) = str2num(DAT1genotypesforMatlab{i,2});
                    DAT1_nosplit(s) = max(str2num(DAT1genotypesforMatlab{i,2}),1);
                else
                    DAT1_split(s) = NaN;
                    DAT1_nosplit(s) = NaN;
                    dud_temp = [dud_temp,s];
                end
            end
        end
    end
    DAT1_tags = {'0/1 DAT1 10-repeats','2 DAT1 10-repeats'};
end
%%
side_tags = {'Left','Right'};
%% Get rid of duds/include only particular subjects

large_CPP = {};
large_N2c = {};

for s2 = 1:length(subject_folder)
    for s = 1:length(Other_System)
        if strcmp(Other_System{s},subject_folder{s2})
            TCD_index(s) = s2;
        end
    end
    for s = 1:length(Monash_System)
        if strcmp(Monash_System{s},subject_folder{s2})
            Monash_index(s) = s2;
        end
    end
    for s = 1:length(large_CPP)
        if strcmp(large_CPP{s},subject_folder{s2})
            large_CPP_index(s) = s2;
        end
    end
    for s = 1:length(large_N2c)
        if strcmp(large_N2c{s},subject_folder{s2})
            large_N2c_index(s) = s2;
        end
    end
end
%%
duds = [];
single_participants = [];
%%
if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    if DAT1_genotypes
        DAT1_split([duds]) = [];
        DAT1_nosplit([duds]) = [];
    end
    EEG_System([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    if DAT1_genotypes
        DAT1_split = DAT1_split(single_participants);
        DAT1_nosplit = DAT1_nosplit(single_participants);
    end
    EEG_System = EEG_System(single_participants);
end

%% Use Current Source Density transformed erp? 1=yes, 0=no
CSD=0;

%% Define channels, having combined Brain Products and Biosemi data


plot_chans = [1:32];
exclude_chans = [];

tester = zeros(33,1);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','labels','plotchans',plot_chans);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','numbers','plotchans',plot_chans);

%% Triggers

%     'Trigger 1: coherence 90, patch 1, motion dir 90, coh motion onset 1800'
%     'Trigger 2: coherence 90, patch 1, motion dir 270, coh motion onset 1800'
%     'Trigger 3: coherence 90, patch 2, motion dir 90, coh motion onset 1800'
%     'Trigger 4: coherence 90, patch 2, motion dir 270, coh motion onset 1800'
%     'Trigger 5: coherence 90, patch 1, motion dir 90, coh motion onset 2800'
%     'Trigger 6: coherence 90, patch 1, motion dir 270, coh motion onset 2800'
%     'Trigger 7: coherence 90, patch 2, motion dir 90, coh motion onset 2800'
%     'Trigger 8: coherence 90, patch 2, motion dir 270, coh motion onset 2800'
%     'Trigger 9: coherence 90, patch 1, motion dir 90, coh motion onset 3800'
%     'Trigger 10: coherence 90, patch 1, motion dir 270, coh motion onset 3800'
%     'Trigger 11: coherence 90, patch 2, motion dir 90, coh motion onset 3800'
%     'Trigger 12: coherence 90, patch 2, motion dir 270, coh motion onset 3800'

% side,motion,ITI
% targcodes = zeros(2,2,3); %DN:([side [1=left 2=right] x motion [1=up 2=down] x ITI [inter-target-interval, 3 levels)
targcodes(1,1,:) = [101 105 109]; % left patch, up motion
targcodes(1,2,:) = [102 106 110]; % left patch, down motion
targcodes(2,1,:) = [103 107 111]; % right patch, up motion
targcodes(2,2,:) = [104 108 112]; % right patch, down motion

%%
fs=500;
numch=32;
rtlim=[0.300 1.500];
% rtlim=[0.300 1.200];

ch_CPP = [25];

ch_lr{1} = [23];
ch_lr{2} = [27];
ch_rl{1} = [27];
ch_rl{2} = [23];

% stim-locked erps
% ts = -0.500*fs:1.800*fs;
% t = ts*1000/fs;
ts = -0.700*fs:1.800*fs;
t = ts*1000/fs;

% resp-locked erps
trs = [-.700*fs:fs*.100];
tr = trs*1000/fs;

BL_erp = [-100,0];

% zscore threshold
z_thresh = 3;

%% Start loop
for s=1:length(allsubj)
    
    if ismember(subject_folder{s},Monash_System)
        EEG_System(s) = 1;
    elseif ismember(subject_folder{s},Other_System)
        EEG_System(s) = 2;
    else
        keyboard
    end
    
    if ismember(subject_folder{s},Other_System)
        EEG_System(s) = 1;
    end
    if ismember(subject_folder{s},Monash_System)
        EEG_System(s) = 2;
    end
    pause(1)
    load([path_temp subject_folder{s} '\' allsubj{s} 'big_dots_erp'],'erp_LPF_8Hz','erp_LPF_35Hz','erp_LPF_35Hz_CSD','allRT','allrespLR','allTrig','allblock_count',...
        'BL_resp_artrej','ET_BL_resp_artrej');
    
    if strcmp(subject_folder{s},'331M_CL') % really odd tiny artifact meant this trial was messing with CSD!
        allRT(53) = 0; allrespLR(53) = 0; allTrig(53) = 0;
    end
    
    if CSD
        erp=double(erp_LPF_35Hz_CSD);
    else
        erp=double(erp_LPF_35Hz);
    end
    
    % Baseline erp
    baseline_erp = mean(erp(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
    erp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp
    
    disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(find(allTrig)))])
    
    erpr = zeros(size(erp,1),length(tr),size(erp,3));
    
    validrlock = zeros(1,length(allRT)); % length of RTs.
    for n=1:length(allRT);
        [blah,RTsamp] = min(abs(t*fs/1000-allRT(n))); % get the sample point of the RT.
        if RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t) & allRT(n)>0 % is the RT larger than 1st stim RT point, smaller than last RT point.
            erpr(:,:,n) = erp(:,RTsamp+trs,n);
            validrlock(n)=1;
        end
    end
    
    
    for side = 1:2
        for iti = 1:3
            for motion = 1:2
                % calcs the indices of the triggers for each
                % appropriate trial type.
                %             conds{s,iti,side} = find(allTrig==targcodes(side,motion, iti) & allrespLR==1 & ...
                %                 allRT>rtlim(1)*fs & allRT<rtlim(2)*fs);
                conds{s,side,motion, iti} = find(allTrig==targcodes(side,motion, iti) & allrespLR==1 & ...
                    allRT>rtlim(1)*fs & allRT<rtlim(2)*fs & BL_resp_artrej==1  & validrlock);
                
                RTs{s,side,motion, iti} = allRT([conds{s,side,motion, iti}])*1000/fs;
                RTs_log{s,side,motion, iti} = log(allRT([conds{s,side,motion, iti}])*1000/fs);
                RT_zs{s,side,motion, iti} = zscore([RTs_log{s,side,motion, iti}]);
                RT_factors(s,side,motion, iti) = mean([RTs{s,side,motion, iti}]); % can't include AR_08_04_14 & MH_14_04_14 because of mistake
                
                hit{s,side,motion, iti} = find(allTrig==targcodes(side,motion, iti) & allrespLR==1 & ...
                    allRT>rtlim(1)*fs & allRT<rtlim(2)*fs);
                miss{s,side,motion, iti} = find(allTrig==targcodes(side,motion, iti) & allrespLR==3);
                
            end
        end
        RT_all(s,side) = mean([RTs{s,:,side,:}]);
        temp = [RTs{s,:,side, :}]; tempz = [RT_zs{s,:,side, :}];
        temp = temp(find(tempz<z_thresh));
        RT_all_zs(s,side) = mean(temp);
        RT_median_all(s,side) = median([RTs{s,:,side, :}]);
        RT_log_all(s,side) = mean([RTs_log{s,:,side, :}]);
        
        hit_all(s,side) = length([hit{s,:,side, :}]);
        miss_all(s,side) = length([miss{s,:,side, :}]);
        hit_rate(s,side) = 100*hit_all(s,side)/(miss_all(s,side)+hit_all(s,side));
    end
    RT_index(s) = (RT_all(s,1)-RT_all(s,2))/((RT_all(s,1)+RT_all(s,2))/2);
    RT_median_index(s) = (RT_median_all(s,1)-RT_median_all(s,2))/((RT_median_all(s,1)+RT_median_all(s,2))/2);
    RT_log_index(s) = (RT_log_all(s,1)-RT_log_all(s,2))/((RT_log_all(s,1)+RT_log_all(s,2))/2);
    RT_index_zs(s) = (RT_all_zs(s,1)-RT_all_zs(s,2))/((RT_all_zs(s,1)+RT_all_zs(s,2))/2);
    
    disp(['Subject ',allsubj{s},' Total Valid Trials: ',num2str(length([conds{s,:,:}])), ...
        ' = ',num2str(round(100*(length([conds{s,:,:}]))/(16*18))),'%'])
    
    for side = 1:2
        ERP_side(s,:,:,side) = squeeze(mean(erp(1:numch,:,[conds{s,:,side, :}]),3));
        ERPr_side(s,:,:,side) = squeeze(mean(erpr(1:numch,:,[conds{s,:,side, :}]),3));
        CPP_side(s,:,side) = squeeze(mean(mean(erp(ch_CPP,:,[conds{s,:,side, :}]),1),3));
        CPPr_side(s,:,side) = squeeze(mean(mean(erpr(ch_CPP,:,[conds{s,:,side, :}]),1),3));
        N2c_side(s,:,side) = squeeze(mean(mean(erp(ch_lr{side},:,[conds{s,:,side, :}]),1),3));
        N2i_side(s,:,side) = squeeze(mean(mean(erp(ch_rl{side},:,[conds{s,:,side, :}]),1),3));
        %% Code adapted from Ger's Current Biology cpp code to pull out CPP onset latency:
        % Define CPP onset search window, from 0 to 1000ms
        CPP_search_t  = [0,1000];
        % Same window in samples
        CPP_search_ts  = [find(t==CPP_search_t(1)),find(t==CPP_search_t(2))];
        % Size of sliding window. This is in fact 1/4 of the search window in ms.
        % So 25 is 100ms. (25 samples x 2ms either side of a particular sample).
        max_search_window = 25;
        
        
        %         consecutive_windows=10;%Number of consecutive windows that p must be less than .05 for in order to call it a CPP onset
        if any(strcmp(subject_folder(s),{'ND_16_05_14'})) || any(strcmp(subject_folder(s),{'036M_JK'}))
            consecutive_windows=50;%had to make it longer for these participants otherwise it records a false CPP onset
        else
            consecutive_windows=15;%15 works well for everybody else
        end
        %%
        CPP_temp = squeeze(mean(erp(ch_CPP,:,[conds{s,:,side, :}]),1)); % time x trial
        CPPs(:,side) = squeeze(mean(CPP_temp(:,:),2)); % average across trial for plot later on, not used to find onsets.
        % constrain the search window according to parameters above.
        CPP_temp = squeeze(CPP_temp(find(t>=CPP_search_t(1) & t<=CPP_search_t(2)),:));
        prestim_temp = find(t<CPP_search_t(1)); % so we can add it on after getting max peak.
        
        % we want sliding windows for each trial, create smoothed waveform.
        clear win_mean win_mean_inds tstats ps
        for trial = 1:size(CPP_temp,2)
            counter = 1;
            for j = max_search_window:2:size(CPP_temp,1)-max_search_window
                win_mean(counter,trial) = mean(CPP_temp([j-max_search_window+1:j+max_search_window-1],trial));
                win_mean_inds(counter) = j;
                counter = counter+1;
            end
        end
        
        % do t-test to zero across the smoothed trials.
        for tt = 1:size(win_mean,1)
            
            if strcmp( subject_folder(s),'AD48C') %This participant has strainge CPP baseline, so do CPP onset t-test against -1.5 instead of against 0
                [~,P,~,STATS] = ttest(win_mean(tt,:),-1.5);
            else
                [~,P,~,STATS] = ttest(win_mean(tt,:));
            end
            tstats(tt) = STATS.tstat;
            ps(tt) = P;
        end
        
        % when does the ttest cross 0.05? If at all?
        %         onsetp05 = find(ps<0.05 & tstats>0,1,'first');
        
        %DN: added this in to explicitly make sure the "consecutive_windows" number of following p-values from onset are also lower than 0.05.
        clear allp05
        allp05= find(ps<0.05 & tstats>0);
        onsetp05=[];
        for i = 1:length(allp05)
            if  (i+consecutive_windows-1)<=length(allp05)
                if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least 10 consecutive windows where p<.05
                    onsetp05=allp05(i);
                    break
                end
            else
                onsetp05=allp05(i);
                break
            end
        end
        
        
        % get timepoint of min index.
        if ~isempty(onsetp05)
            onset_ind = win_mean_inds(onsetp05);
            CPP_onset_ind = onset_ind + length(prestim_temp); % see above, this needs to be added to get the overall time with respect to t.
            CPP_side_onsets(s,side) = t(CPP_onset_ind);
        else % onsetp05 is empty, no significant CPP.
            disp([allsubj{s},': bugger']) %AD48C has no CPP onset
            CPP_side_onsets(s,side) = 0;
        end
        
        % %     plot the smoothed waveforms, the corresponding t-tests and p-values.
        % %     Make sure the 10 (DN:30) following p-values from onset are also lower than
        % %     0.05.
        
        %         if side==1
        %             figure
        %             subplot(3,2,1)
        %             plot(win_mean_inds,mean(win_mean,2))
        %             title(subject_folder{s})
        %             subplot(3,2,3)
        %             plot(win_mean_inds,tstats)
        %             subplot(3,2,5)
        %             plot(win_mean_inds,ps), hold on
        %             line(xlim,[0.05,0.05],'Color','k','LineWidth',1);
        %             if ~isempty(onsetp05)
        %                 line([onset_ind,onset_ind],ylim,'Color','g','LineWidth',1);
        %             else
        %                 line([0,0],ylim,'Color','r','LineWidth',1);
        %             end
        %         else
        %             subplot(3,2,2)
        %             plot(win_mean_inds,mean(win_mean,2))
        %             title(subject_folder{s})
        %             subplot(3,2,4)
        %             plot(win_mean_inds,tstats)
        %             subplot(3,2,6)
        %             plot(win_mean_inds,ps), hold on
        %             line(xlim,[0.05,0.05],'Color','k','LineWidth',1);
        %             if ~isempty(onsetp05)
        %                 line([onset_ind,onset_ind],ylim,'Color','g','LineWidth',1);
        %             else
        %                 line([0,0],ylim,'Color','r','LineWidth',1);
        %             end
        %         end
    end
    %           pause(1)
    %%   plot CPP with onset marked
    %     colors = {'b' 'r' 'g' 'm' 'c'};
    %     figure
    %     for side = 1:2
    %         plot(t,squeeze(CPPs(:,side)),'Color',colors{side, :},'LineWidth',2), hold on
    %         line([mean(CPP_side_onsets(s,side),1),mean(CPP_side_onsets(s,side),1)],ylim,'Color',colors{side, :},'LineWidth',1.5);
    %         line([0,0],ylim,'Color','k','LineWidth',1);
    %         line(xlim,[0,0],'Color','k','LineWidth',1);
    %     end
    %     title(subject_folder{s})
    %     pause(1)
    
    %     erp_mean = (squeeze(mean(mean(mean(ERP_uptarg(s,:,:,:,:,:),4),5),6))+squeeze(mean(mean(mean(ERP_righttarg(s,:,:,:,:,:),4),5),6)))/2;
    %     figure
    %     plottopo(erp_mean(:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
    %         min(min(min(erp_mean(plot_chans,:))))  max(max(max(erp_mean(plot_chans,:))))], ...
    %         'title',[allsubj{s}],'legend',{'90 deg (up)','60 deg','30 deg','0 deg (right)'},'showleg','off','ydir',1)
    %     pause(1)
end

%% Quick check of RT index
[~,p,~,stats] = ttest(RT_index);
disp(['RT index to zero: t = ' num2str(stats.tstat) ', p = ' num2str(p)])

if DAT1_genotypes
[~,p,~,stats] = ttest2(RT_index(find(DAT1_nosplit==1)),RT_index(find(DAT1_nosplit==2)));
disp(['RT index x DAT1: t = ' num2str(stats.tstat) ', p = ' num2str(p)])
end

figure
plot(zscore(RT_index))
%% Grand average ERP
% chan x time x side
ERP_group_side = squeeze(mean(ERP_side(:,:,:,:),1));

figure
plottopo(ERP_group_side(:,:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
    min(min(min(ERP_group_side(plot_chans,:,:))))  max(max(max(ERP_group_side(plot_chans,:,:))))], ...
    'title',['ERP left vs right targets'],'legend',side_tags,'showleg','on','ydir',1)

clear time_windows
time_windows(1,:) = [100:50:500];
time_windows(2,:) = time_windows(1,:)+50;
for side = 1:2
    figure
    for t_temp = 1:size(time_windows,2)
        plot_mean = squeeze(mean(mean(ERP_group_side(:,find(t>time_windows(1,t_temp) & t<time_windows(2,t_temp)),side),2),3));
        subplot(3,3,t_temp)
        topoplot(plot_mean,chanlocs,'maplimits', ...
            [min(min(min(ERP_group_side(:,find(t>time_windows(1,1) & t<time_windows(2,end)),:))))...
            max(max(max(ERP_group_side(:,find(t>time_windows(1,1) & t<time_windows(2,end)),:))))], ...
            'electrodes','off','plotchans',plot_chans);
        title([side_tags{side},' targets: ',num2str(time_windows(1,t_temp)),' ms to ',num2str(time_windows(2,t_temp)),' ms']);
        colorbar
    end
end


%% Make some scalp plots for R:
%CPP scalp topo
t1 = 450; t2 = 600;
plot_mean = squeeze(mean(mean(ERP_group_side(:,find(t>=t1 & t<t2),:),2),3));
figure
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean) max(plot_mean)], ...
    'electrodes','off','plotchans',plot_chans);

%N2 scalp topo
t1 = 200; t2 = 300;
for side = 1:2
    plot_mean = squeeze(mean(mean(ERP_group_side(:,find(t>=t1 & t<t2),side),2),3));
    figure
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean) max(plot_mean)], ...
        'electrodes','off','plotchans',plot_chans);
end


%% Grand average ERPr
% chan x time x side
ERPr_group_side = squeeze(mean(ERPr_side(:,:,:,:),1));

figure
plottopo(ERPr_group_side(:,:,:),'chanlocs',chanlocs,'limits',[tr(1) tr(end) ...
    min(min(min(ERPr_group_side(plot_chans,:,:))))  max(max(max(ERPr_group_side(plot_chans,:,:))))], ...
    'title',['ERPr left vs right targets'],'legend',side_tags,'showleg','on','ydir',1)

clear time_windows
time_windows(1,:) = [-350:50:50];
time_windows(2,:) = time_windows(1,:)+50;
for side = 1:2
    figure
    for t_temp = 1:size(time_windows,2)
        plot_mean = squeeze(mean(mean(ERPr_group_side(:,find(tr>time_windows(1,t_temp) & tr<time_windows(2,t_temp)),side),2),3));
        subplot(3,3,t_temp)
        topoplot(plot_mean,chanlocs,'maplimits', ...
            [min(min(min(ERPr_group_side(:,find(tr>time_windows(1,1) & tr<time_windows(2,end)),:))))...
            max(max(max(ERPr_group_side(:,find(tr>time_windows(1,1) & tr<time_windows(2,end)),:))))], ...
            'electrodes','off','plotchans',plot_chans);
        title([side_tags{side},' targets (resp-locked): ',num2str(time_windows(1,t_temp)),' ms to ',num2str(time_windows(2,t_temp)),' ms']);
        colorbar
    end
end

% t1 = -100; t2 = 0;
% plot_mean = squeeze(mean(mean(ERPr_group_side(:,find(tr>=t1 & tr<t2),:),2),3));
% figure
% topoplot(plot_mean,chanlocs,'maplimits', ...
%     [min(plot_mean) max(plot_mean)], ...
%     'electrodes','numbers','plotchans',plot_chans);
%% Plot CPP x target side
CPP = squeeze(mean(CPP_side,1)); % time x side
clear h
figure
for side = 1:2
    h(side) = plot(t,squeeze(CPP(:,side)),'LineWidth',3,'LineStyle','-');hold on
end
set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('CPP: Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,side_tags, ...
    'FontSize',16,'Location','NorthWest');

% figure
% for s = 1:length(allsubj)
%     CPP = squeeze(mean(CPP_side(s,:,:),3));
%     if max(CPP)>60
%         disp(subject_folder{s})
%     end
%     plot(t,CPP,'LineWidth',1,'LineStyle','-');hold on
% end
% set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
% ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
% xlabel('Time (ms)','FontName','Arial','FontSize',16)
% title('CPP: Left vs Right Hemifield Targets')
% line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
% line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');

%% Plot CPPr x target side
CPPr = squeeze(mean(CPPr_side,1)); % time x side
clear h
figure
for side = 1:2
    h(side) = plot(tr,squeeze(CPPr(:,side)),'LineWidth',3,'LineStyle','-');hold on
end

set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('CPP (resp-locked): Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,side_tags, ...
    'FontSize',16,'Location','NorthWest');

%% Plot N2c x target side
N2c = squeeze(mean(N2c_side,1)); % time x side
clear h
figure
for side = 1:2
    h(side) = plot(t,squeeze(N2c(:,side)),'LineWidth',3,'LineStyle','-');hold on
end
set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('N2c: Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,side_tags, ...
    'FontSize',16,'Location','NorthWest');

% figure
% for s = 1:length(allsubj)
%     N2c = squeeze(mean(N2c_side(s,:,:),3));
%     if min(N2c(find(t<500)))<-5
%         disp(subject_folder{s})
%     else
% %         plot(t,N2c,'LineWidth',1,'LineStyle','-');hold on
%     end
%     plot(t,N2c,'LineWidth',1,'LineStyle','-');hold on
% end
% set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
% ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
% xlabel('Time (ms)','FontName','Arial','FontSize',16)
% title('CPP: Left vs Right Hemifield Targets')
% line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
% line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');

%% Plot N2i x target side
N2i = squeeze(mean(N2i_side,1)); % time x side
clear h
figure
for side = 1:2
    h(side) = plot(t,squeeze(N2i(:,side)),'LineWidth',3,'LineStyle','-');hold on
end
set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('N2i: Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,side_tags, ...
    'FontSize',16,'Location','NorthWest');

% return
%% CPP x DAT1
if DAT1_genotypes
colors = {'b','r'};
linestyles = {'-','--'};
cc=1;
clear h
figure
for dat = 1:2
    CPP = squeeze(mean(CPP_side(find(DAT1_nosplit==dat),:,:),1)); % time x side
    for side = 1:2
        h(cc) = plot(t,squeeze(CPP(:,side)),'Color',colors{side},'LineWidth',2,'LineStyle',linestyles{dat});hold on
        cc=cc+1;
    end
end
set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('CPP: Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,{'Left/DAT1 no repeat','Right/DAT1 no repeat','Left/DAT1 repeat','Right/DAT1 repeat'},...
    'FontSize',16,'Location','NorthWest');

t1 = 270; t2 = 470;
stat_temp = squeeze(mean(CPP_side(:,find(t>=t1 & t<=t2),:),2)); % subj x side
stat_temp = [stat_temp,DAT1_nosplit'];

%% CPPr x DAT1
colors = {'b','r'};
linestyles = {'-','--'};
cc=1;
clear h
figure
for dat = 1:2
    CPPr = squeeze(mean(CPPr_side(find(DAT1_nosplit==dat),:,:),1)); % time x side
    for side = 1:2
        h(cc) = plot(tr,squeeze(CPPr(:,side)),'Color',colors{side},'LineWidth',2,'LineStyle',linestyles{dat});hold on
        cc=cc+1;
    end
end
set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('CPP (resp-locked): Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,{'Left/DAT1 no repeat','Right/DAT1 no repeat','Left/DAT1 repeat','Right/DAT1 repeat'},...
    'FontSize',16,'Location','NorthWest');

%% N2c x DAT1
colors = {'b','r'};
linestyles = {'-','--'};
cc=1;
clear h
figure
for dat = 1:2
    N2c = squeeze(mean(N2c_side(find(DAT1_nosplit==dat),:,:),1)); % time x side
    for side = 1:2
        h(cc) = plot(t,squeeze(N2c(:,side)),'Color',colors{side},'LineWidth',2,'LineStyle',linestyles{dat});hold on
        cc=cc+1;
    end
end
set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('N2c: Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,{'Left/DAT1 no repeat','Right/DAT1 no repeat','Left/DAT1 repeat','Right/DAT1 repeat'},...
    'FontSize',16,'Location','NorthWest');

t1 = 250; t2 = 290;
stat_temp = squeeze(mean(N2c_side(:,find(t>=t1 & t<=t2),:),2)); % subj x side
stat_temp = [stat_temp,DAT1_nosplit'];

%% N2i x DAT1
colors = {'b','r'};
linestyles = {'-','--'};
cc=1;
clear h
figure
for dat = 1:2
    N2i = squeeze(mean(N2i_side(find(DAT1_nosplit==dat),:,:),1)); % time x side
    for side = 1:2
        h(cc) = plot(t,squeeze(N2i(:,side)),'Color',colors{side},'LineWidth',2,'LineStyle',linestyles{dat});hold on
        cc=cc+1;
    end
end
set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title('N2i: Left vs Right Hemifield Targets')
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,{'Left/DAT1 no repeat','Right/DAT1 no repeat','Left/DAT1 repeat','Right/DAT1 repeat'},...
    'FontSize',16,'Location','NorthWest');

end
%% Extract Response Locked CPP slope"
%CPP build-up defined as the slope of a straight line fitted to the
%response-locked waveform at during "slope_timeframe_index" defined for
%each participant here:
clear slope_timeframe_index
for s=1:length(allsubj)
    slope_timeframe_index(s,2)=find(mean(CPPr_side(s,:,:),3)==max(mean(CPPr_side(s,find(tr<0),:),3)));%max amplitude index
end
slope_timeframe_index(:,1)=slope_timeframe_index(:,2)-50;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
%Now find and save CPPr slope
for s=1:length(allsubj)
    for side = 1:2
        coef = polyfit(tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(CPPr_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
        CPPr_slope(s,side)=coef(1);
    end
end

%%%Plot each individual participant's CPPr_slope with time-window varying
%%%per participant
for s=1:length(allsubj)
    clear h
    figure
    for side = 1:2
        h(side) = plot(tr,CPPr_side(s,:,side),'LineWidth',3,'LineStyle','-');hold on
        coef = polyfit(tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(CPPr_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
        CPP_slope(s,side)=coef(1);
        r = coef(1) .* tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
        plot(tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)), r,'Linewidth',2, 'LineStyle', ':');
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder{s}, ' CPP (resp-locked) by Hemifield'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    pause(1)
end



%% Extract N2c and N2i latency :
for s = 1:size(allsubj,2)
    for side=1:2
        %         to use each participant's average N2c/i to get their peak latency index:
        N2c=N2c_side(s, :, side);
        N2i=N2i_side(s, :, side);
        N2c_peak_amp_index_t(s,side)=t(N2c==min(N2c(find(t==150):find(t==400))));%Find max peak latency for N2c in ms
        N2i_peak_amp_index_t(s,side)=t(N2i==min(N2i(find(t==200):find(t==450))));%Find max peak latency for N2i in ms
    end
end
%N2 Latency:
N2cN2i_latency_ByTargetSide = [N2c_peak_amp_index_t,N2i_peak_amp_index_t]; %(LeftTargetN2c_latency, RightTargetN2c_latency, LeftTargetN2i_latency, RightTargetN2i_latency)

% %Plot N2c and N2i per subject showing peak amplitude
% for s = 1:size(allsubj,2)
% clear h
% figure
% for side = 1:2
%     h(side) = plot(t,squeeze(N2i_side(s,:,side)),'LineWidth',3,'LineStyle',':', 'Color',colors{side});hold on
%     h(side) = plot(t,squeeze(N2c_side(s,:,side)),'LineWidth',3,'LineStyle','-', 'Color',colors{side});hold on
%     line([avN2i_peak_amp_index_t(s,side),avN2i_peak_amp_index_t(s,side)],ylim,'Color',colors{side}, 'Linewidth',2,'LineStyle',':');
%     line([avN2c_peak_amp_index_t(s,side),avN2c_peak_amp_index_t(s,side)],ylim,'Color',colors{side}, 'Linewidth',2,'LineStyle','-');
% end
%
% set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
% ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
% xlabel('Time (ms)','FontName','Arial','FontSize',16)
% title(['Subj: ',num2str(s),' N2c & N2i(:) by Hemifield'])
% line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
% line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
% legend(h,side_tags, 'FontSize',16,'Location','NorthWest');
% pause(1)
% end

%% Extract N2c and N2i Amplitude :
window=25; %this is the time (in samples) each side of the peak latency (so it's 50ms each side of peak latency - so a 100ms window)
N2c = squeeze(mean(mean(N2c_side,1),3)); % time
N2i = squeeze(mean(mean(N2i_side,1),3)); % time
N2c_peak_amp_index=find(N2c==min(N2c(find(t==150):find(t==450))));%Find Left target max peak latency for N2c
N2i_peak_amp_index=find(N2i==min(N2i(find(t==150):find(t==450))));%Find Left target max peak latency for N2i

for s = 1:size(allsubj,2)
    for side=1:2
        max_peak_N2c(s,side)=squeeze(mean(N2c_side(s,N2c_peak_amp_index-window:N2c_peak_amp_index+window, side),2));
        max_peak_N2i(s,side)=squeeze(mean(N2i_side(s,N2i_peak_amp_index-window:N2i_peak_amp_index+window, side),2));
    end
end
N2cN2i_amp_ByTargetSide_ParticipantLevel = [max_peak_N2c,max_peak_N2i]; %(LeftTargetN2c, RightTargetN2c, LeftTargetN2i, RightTargetN2i)


%% Make participant level matrix for export into SPSS or R
participant_level(:,1:2)=max_peak_N2c; %N2c amplitude (LeftTarget, RightTarget)
participant_level(:,3:4)=max_peak_N2i; %N2i amplitude (LeftTarget, RightTarget)
participant_level(:,5:6)=N2c_peak_amp_index_t; %N2c latency (LeftTarget, RightTarget)
participant_level(:,7:8)=N2i_peak_amp_index_t; %N2i latency (LeftTarget, RightTarget)
participant_level(:,9:10)=CPP_side_onsets; %CPP onset (LeftTarget, RightTarget)
participant_level(:,11:12)=CPPr_slope; %response locked CPP slope (LeftTarget, RightTarget)
participant_level(:,13)=EEG_System;%1=Monash_System; 2=Other_System
% open participant_level

csvwrite (['participant_level_matrix.csv'],participant_level)

subject_folder=subject_folder';
cell2csv ('IDs.csv',subject_folder)