
function rm=returnmap(spikes,plotflag)
if size(spikes,2)>1, spikes=spikes';end
isis=diff(spikes); isis=[spikes(1);isis];
rm.d1=isis(1:end-1);%preceding ISIs
rm.d2=isis(2:end); %followoing isis
rm.spikes=isis(1:end-1);%preceding ISIs
rm.r=corrcoef([rm.d1 rm.d2]);

if plotflag
    rm.fig=figure;
    plot(log10(rm.d1),log10(rm.d2),'.')
    xlabel preceding
    ylabel following
end

