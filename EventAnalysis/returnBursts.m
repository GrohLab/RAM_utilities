
function [Bursts firstSpikes C BurstLength]=returnBursts(sp,isiCutoff)
%[Bursts firstSpikes C]=returnBursts(sp,isiCutoff)
%This function takes a spike train, sp, and an interspike interval
%criterion, isiCutoff, and returns a structure of events, Bursts; each
%entry represents a burst of spikes sorted according to isiCutoff.
%firstSpikes gives a vector of the first spike in each event.
%C gives a vector of each event's size;
%BurstLength--time from first spike to last spike
%numel(C)=numel(firstSpikes)=numel(Bursts);
BurstLength=[];
if numel(sp)==0
    Bursts=[];
    firstSpikes=[];
    C=[];
end
% what if there is just one spike?
if numel(sp)==1
    Bursts={sp}
    firstSpikes=sp;
    C=1;
else
    %if there are multiple spikes, then sort them;
    dim=size(sp); if dim(2)>dim(1),sp=sp';end  %consistent column of spike times.
    Bursts={};
    %what if there are no spikes?
    if isempty(sp)
        Bursts={};
    elseif numel(sp)==1   %this is redundant, but left for now.
        Bursts{1}=sp;
    else
        isis=[sp(1); diff(sp)];  %find interspike intervals
        bursts=(find(isis>isiCutoff)); %indices of burst starts;
        bursts=[1;bursts];bursts=unique(bursts);% what if first even is too close to beginning? count it as event.
        firstSpikes=sp(bursts); %find first spikes;
        C=nan(size(bursts)); %preallocate
        for i=1:numel(bursts)
            % for each burst starting spike, find spikes between it and the
            % the next burst starting spike.
            if i<numel(bursts)
                i1=bursts(i);
                i2=bursts(i+1)-1;
                Bursts{i}=[sp(i1:i2)];
            else
                i1=bursts(i);
                Bursts{i}=[sp(i1:end)];
            end
            C(i)=numel(Bursts{i});  %how many spikes per burst?
        end
    end
end
if ~isempty(Bursts)
    for i=1:numel(Bursts)
        BurstLength(i)=sum(diff(Bursts{i}));
    end
end

