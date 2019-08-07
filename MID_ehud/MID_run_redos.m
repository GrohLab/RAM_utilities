
cd C:\Users\rebecca\Dropbox\SynchronyFall2014
clear all
%load all 23 single cells
%script
GetSinglePOm_cells
DS=20; %downsampling for eeg preprocessing (to 1000 Hz)

%%

%  run events only with new parameters, parallelized for speed
pack
MID={};
Ibest={};
Inf={};
stas={}
FirstSpikes={};
cd C:\Users\rebecca\Dropbox\synchronyFall2014\MID_ehud
timeBefore=.4
timeAfter=0;
ds=1000;
%grab variables for parallel run
fs=500
for I=1:23
    spikes=Data{I}.spikes;
    %extract bursts
    [Bursts firstSpikes Cs{I}]=returnBursts(spikes,20*20);
    spikes=firstSpikes/20000; %first spikes of events, in seconds
    eeg=Data{I}.eeg_clean; eeg=zscore(eeg);
    EEGs{I}=eeg(2:2:end); %downsample to 500 Hz
    FirstSpikes{I}=spikes;
end


EEGs=EEGs(sort([22 21 18 17 16 12 11 10]))
FirstSpikes=FirstSpikes(sort([22 21 18 17 16 12 11 10]))




%%
MID={};
Ibest={};
Inf={};
stas={}
%% redos
parfor I=1:8
    display(I)
    [t,MID{I},Ibest{I},Inf{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,0);
   % [t,MID_sta{I},Ibest_sta{I},Inf_sta{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,1);
    %     save POM_mid_events_t6 MID Ibest Inf stas MID_sta Ibest_sta Inf_sta Cs FirstSpikes t
end
save POM_mid_events_t10_redo MID Ibest Inf stas   



%%
MIDredo={};

for i=1:8
    MIDredo{indices(i)}=MID{i};
end


%%
Data=load('POM_mid_events_t6')

MID=Data.MID;
for i=indices
    plot(Data.MID{i});
    hold on
    try
    plot(MIDredo{i})
    catch
    end
    pause
    hold off
    MID{i}=MIDredo{i};
end
    


%% new redo (start with best guess from last time)

% 12 16 21 22  goes to redos   3 4 7 8
load POM_mid_events_t10_redo MID

guess=MID([3 4 7 8])

%%
%  run events only with new parameters, parallelized for speed
pack
MID={};
Ibest={};
Inf={};
stas={}
FirstSpikes={};
cd C:\Users\rebecca\Dropbox\synchronyFall2014\MID_ehud
timeBefore=.4
timeAfter=0;
ds=1000;
%grab variables for parallel run
fs=500
for I=1:23
    spikes=Data{I}.spikes;
    %extract bursts
    [Bursts firstSpikes Cs{I}]=returnBursts(spikes,20*20);
    spikes=firstSpikes/20000; %first spikes of events, in seconds
    eeg=Data{I}.eeg_clean; eeg=zscore(eeg);
    EEGs{I}=eeg(2:2:end); %downsample to 500 Hz
    FirstSpikes{I}=spikes;
end
EEGs=EEGs(sort([12 16 21 22]))
FirstSpikes=FirstSpikes(sort([12 16 21 22]))



%%
MID={};
Ibest={};
Inf={};
stas={}
%% redos
parfor I=1:4
    display(I)
    [t,MID{I},Ibest{I},Inf{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,guess{I});
   % [t,MID_sta{I},Ibest_sta{I},Inf_sta{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,1);
    %     save POM_mid_events_t6 MID Ibest Inf stas MID_sta Ibest_sta Inf_sta Cs FirstSpikes t
end
save POM_mid_events_t10_redo2 MID Ibest Inf stas guess



%% redos
parfor I=1:2
    display(I)
    [t,MID{I},Ibest{I},Inf{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,1);
   % [t,MID_sta{I},Ibest_sta{I},Inf_sta{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,1);
    %     save POM_mid_events_t6 MID Ibest Inf stas MID_sta Ibest_sta Inf_sta Cs FirstSpikes t
end
save POM_mid_events_t10_redo3_sta MID Ibest Inf stas guess



%%

figure
for i=1:numel(MID)
    plot(MID{i})
    hold on
    plot(guess{i})
    hold off
    pause
end

counter=0
for i=([12 16 22 21])
    counter=counter+1
    MID{i}=finalMIDredo{counter};
end
