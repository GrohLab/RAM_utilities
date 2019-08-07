
function spikes=getspiketimes2(v,thresh,spikewidth)
        
spikes=find(v>=thresh);

spikes2=round(spikes/spikewidth)*spikewidth;

windows=unique(spikes2);
spiketimes=[];
for i=1:numel(windows)
    time=find(spikes2==windows(i));
    spike=spikes(time);
    spiketime=spike(find(v(spike)==max(v(spike))));
    spiketime=spiketime(1);
    if( spiketime>1 &&  v(spiketime-1)<v(spiketime) && spiketime<numel(v) && v(spiketime)>v(spiketime+1))
        spiketimes=[spiketimes; spiketime];
    end
end

spikes=spiketimes;

 