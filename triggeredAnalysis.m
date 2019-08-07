function [sp h bins uptrig_EEG f]=triggeredAnalysis(spikes,ppms,triggers,binsize,timeBefore, timeAfter,str, EEG, plotit);

f=0;
if size(spikes,2)>size(spikes,1),spikes=spikes';end
sp={};
for i=1:numel(triggers)
    spnorm=spikes-triggers(i);
    sp{i}=spnorm(find(spnorm>=-timeBefore & spnorm<=timeAfter));
end

if plotit
    f=figure
    subplot(6,1,1:4)    
    for i=1:numel(sp)
        if ~isempty(sp{i})
            r=raster(sp{i}/ppms,i,'k',.9)
        end
        hold on
        
    end
    try
        set(r,'ShowBaseLine','off');
    catch
    end
    ylabel('trial')
    ylim([0 numel(sp)])
    xlim([-timeBefore/ppms  timeAfter/ppms])
    
    
end
bins=[-timeBefore:binsize:timeAfter]/ppms;
h=hist(cell2mat(sp')/ppms,bins)/numel(sp);


if plotit
    plot([0 0], [0 numel(sp)+1],'r:','linewidth',3)
    title(str)
    subplot(6,1,5)
    
    
    bar(bins,h,1);
    xlim([-timeBefore/ppms  timeAfter/ppms])
    hold on
    plot([0 0], get(gca,'ylim'),'r:','linewidth',3)
    
    ylabel('count/trial')
    
    subplot(6,1,6)
    
end


uptrig_EEG=0;
if ~isempty(EEG)
    uptrig_EEG=TriggeredSegments(EEG,triggers,-timeBefore, timeAfter)';
end

if plotit
    t=[-timeBefore:timeAfter]/ppms;
    p=plot(t,mean(uptrig_EEG),'k','linewidth',2)
    hold on
    axis tight
    ylim(get(gca,'ylim')*1.5);
    plot([0 0], get(gca,'ylim'),'r:','linewidth',3)
    xlabel('ms relative to trigger')
    set(gcf,'position',[597    49   483   636])
end
