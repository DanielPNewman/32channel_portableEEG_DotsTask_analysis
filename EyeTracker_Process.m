        % ET stuff
        load(ET_matfiles{f}) %DN: load the ET mat file
        
        % fix missing trig issue
        EEG_trigs=[]; ET_trigs=[];
        for i = 1:length(EEG.event), EEG_trigs(i) = EEG.event(i).type; end
        for i = 1:length(event), ET_trigs(i) = event(i,2); end
        ET_trigs = ET_trigs(find(ET_trigs>100)); EEG_trigs = EEG_trigs(find(EEG_trigs>100));
        
        if length(ET_trigs)>length(EEG_trigs), last_event = ET_trigs(end-1);
            if ET_trigs(end)==ET_trigs(end-1), event = event(1:end-2,:); save(ET_matfiles{f},'event','-append'); end
        end
        
        plotter = 0;
        %Add an extra 4 rows into the EEG struct - 'TIME'
        %'GAZE_X' 'GAZE_Y' 'AREA'. This will add these as extra channels onto %EEG.data. So the final channel is the pupil area (i.e. diameter):
        EEG_temp = pop_importeyetracker(EEG,ET_matfiles{f},[first_event ...
            last_event],[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},0,1,0,plotter);
        [output_cell,~,~] = command_window_text();
        text = output_cell{length(output_cell)-1};
        numsamp = sscanf(text,'%*s%*s%*s%*s%*s%*f%*s%*s%f');
        if numsamp<30
            beep
            disp([allsubj{s},', block ',num2str(f),': ET sync issue'])
            figure, plot(ET_trigs), hold on, plot(EEG_trigs)
            keyboard
        end
    
    ET_data = EEG_temp.data([nchan+1:nchan+4],:);
%     scres = 1280 x 1024: 640,512 is middle
%     temp = ET_data([2,3],[stimes(1):end]); ET_data_mean(1) = mean(temp(1,find(temp(1,:)>0)),2); ET_data_mean(2) = mean(temp(2,find(temp(2,:)>0)),2);
%     ET_data(2,find(ET_data(2,:)>0)) = ET_data(2,find(ET_data(2,:)>0))-ET_data_mean(1)+640;
%     ET_data(3,find(ET_data(3,:)>0)) = ET_data(3,find(ET_data(3,:)>0))-ET_data_mean(2)+512;
    clear EEG_temp
    pause(1)
