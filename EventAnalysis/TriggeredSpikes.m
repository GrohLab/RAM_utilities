function [TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(Spikes,triggers,timeBefore, timeAfter)

N_Sp=numel(Spikes);n_trig=numel(triggers); %number of units, number of triggers
TriggeredSpikeTimes=cell(1,N_Sp); %for each neuron or unit
Yvals={};Xvals={}; %make structures to collect data for plotting and psth creation;
shift=0;
for I= 1:N_Sp  %per neuron
    spikes=double(Spikes{I});%spike times
    if size(spikes,2)>size(spikes,1),spikes=spikes';end
    
    TriggeredSpikeTimes{I}=cell(1,n_trig);
    ys={};
    for i=1:n_trig   %for each trigger
        spnorm=spikes-triggers(i); %adjust to be relative to trigger
        TriggeredSpikeTimes{I}{i}=spnorm(find(spnorm>=-timeBefore & spnorm<=timeAfter)); %only spikes in window
        
        %plotting values
        if ~isempty(TriggeredSpikeTimes{I}{i})
            ys{i}=ones(size(TriggeredSpikeTimes{I}{i}))+shift;
        else
            ys{i}=[];
        end
        shift=shift+1;%add for each trial, per neuron
    end
    Xvals{I}=cell2mat(TriggeredSpikeTimes{I}');
    Yvals{I}=cell2mat(ys');
end

