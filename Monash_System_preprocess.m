
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


%% Time

targcodes = [101:112];
nchan = 32; %32. 33 channels if you include the reference channel FCz is added back in at channel 33)

fs = 500; % new sample rate

% ts = -0.500*fs:1.800*fs;
% t = ts*1000/fs;
ts = -0.700*fs:1.800*fs;
t = ts*1000/fs;

BL_time = [-100 0];   % baseline interval in ms
default_response_time = 1.700;

ERP_samps = length(ts);

%% Filters

LPFcutoff_35Hz=35;       % Low Pass Filter cutoff
LPFcutoff_8Hz=8;       % Low Pass Filter cutoff

HPFcutoff=0.1;       % High Pass Filter cutoff

LPF = 1;    % 1 = low-pass filter the data, 0=don't.
HPF = 0;

%% Artifact rejection

ARchans = [1:32];  % just the ones we care about, and RH opposite LH to be symmetric
artifth = 100;
artifchans=[];  % keep track of channels on which the threshold is exceeded, causing trial rejection

chanlocs = readlocs('actiCAP32_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
chanlocs = chanlocs(1:nchan)';

%% Initialise
% # % Make sure all initialised
numtr=0;
allRT=[]; allrespLR=[]; allTrig=[]; allblock_count = [];
erp_LPF_8Hz = []; erp_LPF_35Hz = []; erp_LPF_8Hz_CSD = []; erp_LPF_35Hz_CSD = [];
artifchans_pretarg = []; artifchans_BL_resp = []; artifchans_resp = []; artifchans_1000ms = [];
pretarg_artrej = []; BL_resp_artrej = []; resp_artrej = []; t1000ms_artrej = [];
ET_pretarg_artrej = []; ET_BL_resp_artrej = []; ET_resp_artrej = []; ET_t1000ms_artrej = []; ET_trials=[];

%% Begin loop
disp(subject_folder{s})
for f=1:length(files)
    disp(f)
    EEG = pop_loadbv(paths{f},files{f},[],1:nchan);
    %     loadbvSK_DN % this takes 33 channels (i.e. the reference channel FCz is added back in at cnannel 33)
    EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
    EEG.data = double(EEG.data);
    while 1
        if EEG.event(1).type==1 | EEG.event(1).type==-88
            EEG.event(1) = [];
        else
            break;
        end
    end
    
    load(matfiles{f},'trialCond','par');
    
    %% Eye-tracker?
    if use_el
        fid = fopen(ET_files{f});
        ET_text = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','Headerlines',22,'ReturnOnError',0);
        fclose(fid);
        for i = 1:size(ET_text{1,3},1)
            if strcmp('GAZE_COORDS',ET_text{1,3}(i))
                screen_res(1) = str2num(cell2mat(ET_text{1,6}(i)))+1
                screen_res(2) = str2num(cell2mat(ET_text{1,7}(i)))+1
                continue
            end
        end
        if screen_res(1)==1024, ranger = 76; elseif screen_res(1)==1280, ranger = 98; else disp(screen_res), keyboard, end
        middle = screen_res(1)/2;
    end
    
    %     if par.cohLevels~=50
    %         disp(['Coh discrepancy: Coh = ',num2str(par.cohLevels)])
    %         keyboard
    %     end
    trialCond = trialCond+100;
    
    if use_el
        if ~exist(ET_matfiles{f}, 'file') %DN: if ET matfile NOT been saved
            FixEyelinkMessages %then calculate and save it now
        end
    end
    
    numev = length(EEG.event);
    % Fish out the event triggers and times
    clear trigs stimes
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end
    
    % interpolate bad channels
    if ~isempty(badchans)
        EEG.chanlocs = chanlocs;
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end
    
    if use_el
        % ET stuff
        EyeTracker_Process
    end
    % Check fft for mains frequencies etc
    %     ep = squeeze(EEG.data(31,:)); % time
    %     nfft = size(ep,2);
    %     fftx = abs(fft(ep,[],2))./(nfft/2);
    %     fftx = fftx(:,1:ceil((nfft+1)/2));
    %     freq_temp = (0:ceil((nfft+1)/2)-1)*fs/nfft;
    %     figure, plot(freq_temp,fftx)
    
    % Get rid of mains frequency
    EEG = pop_eegfiltnew(EEG,49,51,[],1,0,0,0);
    EEG_LPF_8Hz = EEG; EEG_LPF_35Hz = EEG; clear EEG;
    
    % First LP Filter
    if LPF
        EEG_LPF_8Hz = pop_eegfiltnew(EEG_LPF_8Hz,0,LPFcutoff_8Hz,[]);
        EEG_LPF_35Hz = pop_eegfiltnew(EEG_LPF_35Hz,0,LPFcutoff_35Hz,[]);
    end
    
    % First HP Filter
    if HPF
        EEG_LPF_8Hz = pop_eegfiltnew(EEG_LPF_8Hz,HPFcutoff,0,[]); % filter to 0.1Hz
        EEG_LPF_35Hz = pop_eegfiltnew(EEG_LPF_35Hz,HPFcutoff,0,[]); % filter to 0.1Hz
        disp('HPF finished')
    end
    
    % average-reference the whole continuous data (safe to do this now after interpolation and filtering):
    EEG_LPF_8Hz.data = EEG_LPF_8Hz.data - repmat(mean(EEG_LPF_8Hz.data(:,:),1),[nchan,1]);
    EEG_LPF_35Hz.data = EEG_LPF_35Hz.data - repmat(mean(EEG_LPF_35Hz.data(:,:),1),[nchan,1]);
    
    targtrigs = [];
    for n=1:length(trigs)
        if any(targcodes(:)==trigs(n))
            targtrigs = [targtrigs n];
        end
    end
    
    if trigs(targtrigs(end))==trialCond(1)
        motion_on = targtrigs(1:end-1); % GL: indices of trigs when motion on. get rid of last trig, it was a repeat
    else
        motion_on = targtrigs;
    end
    
    for n=1:length(motion_on)
        clear ep_LPF_8Hz ep_LPF_35Hz ep_LPF_8Hz_CSD ep_LPF_35Hz_CSD
        locktime = stimes(motion_on(n));
        if motion_on(n)<length(trigs)
            if trigs(motion_on(n)+1)==12
                response_time = stimes(motion_on(n)+1)-locktime; % time in samples from beginning of motion to response.
                response_time = floor(response_time);
                if response_time>default_response_time*fs
                    response_time = default_response_time*fs;
                end
            else
                response_time = default_response_time*fs;
            end
        else
            response_time = default_response_time*fs;
        end
        try
            ep_LPF_8Hz = EEG_LPF_8Hz.data(:,locktime+ts);   % chop out an epoch
            ep_LPF_35Hz = EEG_LPF_35Hz.data(:,locktime+ts);
            if use_el
                ep_ET = ET_data(:,locktime+ts);
            end
        catch
            disp('EEG ended too soon2')
            numtr = numtr+1;
            allTrig(numtr) = 0;
            targMotion(numtr) = 0;
            allblock_count(numtr) = f;
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
            
            erp_LPF_8Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_LPF_35Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_LPF_8Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_LPF_35Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            if use_el
                ET_trials(:,:,numtr) = zeros(4,ERP_samps);
            end
            
            continue;
        end
        
        % new trial
        numtr = numtr+1;
        
        BLamp = mean(ep_LPF_8Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); % record baseline amplitude for each channel, ONLY FOR ART REJECT
        ep_LPF_8Hz_BL = ep_LPF_8Hz - repmat(BLamp,[1,length(t)]);
        
        BLamp = mean(ep_LPF_35Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); % record baseline amplitude for each channel, ONLY FOR ART REJECT
        ep_LPF_35Hz_BL  = ep_LPF_35Hz - repmat(BLamp,[1,length(t)]);
        
        ep_test = [find(ts<=(response_time+0.1*fs))];
        if isempty(ep_test)
            disp('Empty epoch for art rejection')
            keyboard
        end
        
        artifchans_thistrial_pretarg = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts<=0))),[],2)>artifth));
        artifchans_thistrial_BL_resp = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts>=-0.1 & ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial_resp = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial_1000ms = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts>=-0.1 & ts<=1))),[],2)>artifth));
        
        artifchans_pretarg = [artifchans_pretarg artifchans_thistrial_pretarg];
        artifchans_BL_resp = [artifchans_BL_resp artifchans_thistrial_BL_resp];
        artifchans_resp = [artifchans_resp artifchans_thistrial_resp];
        artifchans_1000ms = [artifchans_1000ms artifchans_thistrial_1000ms];
        
        % 0 = reject, 1 = keep
        if length(artifchans_thistrial_pretarg) > 0, pretarg_artrej(numtr) = 0; else pretarg_artrej(numtr) = 1; end
        if length(artifchans_thistrial_BL_resp) > 0, BL_resp_artrej(numtr) = 0; else BL_resp_artrej(numtr) = 1; end
        if length(artifchans_thistrial_resp) > 0, resp_artrej(numtr) = 0; else resp_artrej(numtr) = 1; end
        if length(artifchans_thistrial_1000ms) > 0, t1000ms_artrej(numtr) = 0; else t1000ms_artrej(numtr) = 1; end
        
        % scres = 1024 x 768: 512, 384 is middle. 3 deg is 76 pixels. Nope!
        if use_el
            artif_ET_pretarg = find(ep_ET(2,find(ts<=0))<middle-ranger | ep_ET(2,find(ts<=0))>middle+ranger);
            artif_ET_BL_resp = find(ep_ET(2,find(ts>=-0.1 & ts<=response_time+0.1*fs))<middle-ranger | ep_ET(2,find(ts>=-0.1 & ts<=response_time+0.1*fs))>middle+ranger);
            artif_ET_resp = find(ep_ET(2,find(ts<=(response_time+0.1*fs)))<middle-ranger | ep_ET(2,find(ts<=(response_time+0.1*fs)))>middle+ranger);
            artif_ET_1000ms = find(ep_ET(2,find(ts>=-0.1 & ts<=1))<middle-ranger | ep_ET(2,find(ts>=-0.1 & ts<=1))>middle+ranger);
            
            
            % 0 = reject, 1 = keep
            if length(artif_ET_pretarg) > 0, ET_pretarg_artrej(numtr) = 0; else ET_pretarg_artrej(numtr) = 1; end
            if length(artif_ET_BL_resp) > 0, ET_BL_resp_artrej(numtr) = 0; else ET_BL_resp_artrej(numtr) = 1; end
            if length(artif_ET_resp) > 0, ET_resp_artrej(numtr) = 0; else ET_resp_artrej(numtr) = 1; end
            if length(artif_ET_1000ms) > 0, ET_t1000ms_artrej(numtr) = 0; else ET_t1000ms_artrej(numtr) = 1; end
        end
        
        ep_LPF_8Hz = double(ep_LPF_8Hz);
        ep_LPF_35Hz = double(ep_LPF_35Hz);
        
        ep_LPF_8Hz_CSD = CSD(ep_LPF_8Hz,G_CSD,H_CSD);
        ep_LPF_35Hz_CSD = CSD(ep_LPF_35Hz,G_CSD,H_CSD);
        
        erp_LPF_8Hz(:,:,numtr) = ep_LPF_8Hz;
        erp_LPF_35Hz(:,:,numtr) = ep_LPF_35Hz;
        erp_LPF_8Hz_CSD(:,:,numtr) = ep_LPF_8Hz_CSD;
        erp_LPF_35Hz_CSD(:,:,numtr) = ep_LPF_35Hz_CSD;
        
        if use_el
            ET_trials(:,:,numtr) = ep_ET(:,:);
        end
        
        allTrig(numtr) = trigs(motion_on(n));
        allblock_count(numtr) = f;
        
        try % change this
            if trigs(motion_on(n)+1)==12
                allrespLR(numtr) = 1;
                allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
            elseif trigs(motion_on(n)+1)==13 % they pressed the wrong button
                allrespLR(numtr) = 2;
                allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
            else
                allrespLR(numtr) = 3; % no response, to mark it out from artifact trials.
                allRT(numtr) = 0;
            end
        catch
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
        end
    end
end

rejected_trials = length(find(BL_resp_artrej==0));
figure;
hist(artifchans_BL_resp,[1:nchan]); title([allsubj{s} ': ' num2str(rejected_trials) ' artifacts = ',num2str(round(100*(rejected_trials/length(allRT)))),'%']) % s from runafew
disp([allsubj{s},' number of trials: ',num2str(length(find(allRT)))])

[counts,centers] = hist(artifchans_BL_resp,[1:nchan]);
figure;
topoplot(counts,chanlocs,'plotchans',[1:nchan],'electrodes','numbers');
title(subject_folder{s})
pause(1)

erp_LPF_8Hz = single(erp_LPF_8Hz);
erp_LPF_35Hz = single(erp_LPF_35Hz);
erp_LPF_8Hz_CSD = single(erp_LPF_8Hz_CSD);
erp_LPF_35Hz_CSD = single(erp_LPF_35Hz_CSD);

%% reorganise chanlocs if combining data collected on different systems
% chanlocs_TCD = readlocs('cap64.loc');
% chanlocs_Monash = readlocs('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
%
% counter = 1; counter2 = 1;
% order_TCD=[]; order_Monash=[]; unmatched_elecs_Monash=[];
% for elec_Monash = 1:65
%     yesser = 0;
%     for elec_TCD = 1:64
%         if strcmp(chanlocs_TCD(1,elec_TCD).labels,chanlocs_Monash(1,elec_Monash).labels)
%             order_TCD(counter) = elec_TCD;
%             order_Monash(counter) = elec_Monash;
%             counter = counter+1;
%             yesser = 1;
%         end
%     end
%     if yesser==0
%         unmatched_elecs_Monash(counter2) = elec_Monash;
%         counter2 = counter2+1;
%     end
% end
% erp_TCD = single(zeros(size(erp_LPF_8Hz)));
% erp_TCD(order_TCD,:,:) = erp_LPF_8Hz(order_Monash,:,:);
% erp_LPF_8Hz = erp_TCD(1:nchan,:,:);
% erp_TCD = single(zeros(size(erp_LPF_35Hz)));
% erp_TCD(order_TCD,:,:) = erp_LPF_35Hz(order_Monash,:,:);
% erp_LPF_35Hz = erp_TCD(1:nchan,:,:);
% erp_TCD = single(zeros(size(erp_LPF_8Hz_CSD)));
% erp_TCD(order_TCD,:,:) = erp_LPF_8Hz_CSD(order_Monash,:,:);
% erp_LPF_8Hz_CSD = erp_TCD(1:nchan,:,:);
% erp_TCD = single(zeros(size(erp_LPF_35Hz_CSD)));
% erp_TCD(order_TCD,:,:) = erp_LPF_35Hz_CSD(order_Monash,:,:);
% erp_LPF_35Hz_CSD = erp_TCD(1:nchan,:,:);

%%

% % Baseline erp
% baseline_erp = mean(erp_LPF_35Hz(:,find(t>=BL_time(1) & t<=BL_time(2)),:),2);
% erp_temp = erp_LPF_35Hz-repmat(baseline_erp,[1,size(erp_LPF_35Hz,2),1]); % baseline full erp
%
% erp_temp = squeeze(mean(erp_temp(:,:,find(BL_resp_artrej==1)),3));
% figure
% plottopo(erp_temp(:,:),'chanlocs',chanlocs_TCD,'limits',[t(1) t(end) ...
%     min(min(erp_temp(:,:)))  max(max(erp_temp(:,:)))], ...
%     'title',['ERP'],'ydir',1)


save([path_temp subject_folder{s} '\' allsubj{s} 'big_dots_erp'],'erp_LPF_8Hz','erp_LPF_35Hz','erp_LPF_8Hz_CSD','erp_LPF_35Hz_CSD', ...
    'allRT','allrespLR','allTrig','allblock_count','t','ET_trials', ...
    'artifchans_pretarg','artifchans_BL_resp','artifchans_resp','artifchans_1000ms', ...
    'pretarg_artrej','BL_resp_artrej','resp_artrej','t1000ms_artrej', ...
    'ET_pretarg_artrej','ET_BL_resp_artrej','ET_resp_artrej','ET_t1000ms_artrej')


% % % % % %% Longer epochs for alpha
% % % % %
% % % % % erp_LPF_35Hz_long = erp_LPF_35Hz;
% % % % % erp_LPF_35Hz_CSD_long = erp_LPF_35Hz_CSD;
% % % % % pretarg_artrej_long = pretarg_artrej;
% % % % % resp_artrej_long = resp_artrej;
% % % % % ET_pretarg_artrej_long = ET_pretarg_artrej;
% % % % % ET_resp_artrej_long = ET_resp_artrej;
% % % % % t_long = t;
% % % % %
% % % % % save([path_temp subject_folder{s} '\' allsubj{s} 'big_dots_erp'],'erp_LPF_35Hz_long','erp_LPF_35Hz_CSD_long', ...
% % % % %     'pretarg_artrej_long','resp_artrej_long','ET_pretarg_artrej_long','ET_resp_artrej_long','t_long','-append')

return;