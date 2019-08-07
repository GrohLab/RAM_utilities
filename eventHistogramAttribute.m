function [h bins trig_events trig_attributes C]=eventHistogramAttribute(events,trigs,N,binsize,timeBefore,timeAfter,attributes)


bievents=zeros(N,1);
biAttributes=nan(N,1);
bievents(events)=1;
biAttributes(events)=attributes;

X=TriggeredSegments(bievents,trigs,timeBefore,timeAfter);
Xa=TriggeredSegments(biAttributes,trigs,timeBefore,timeAfter);

trig_events={}; trig_attributes={};


n=size(X,2);
C=zeros(n,1);
for j=1:n
    trig_events{j}=find(X(:,j))+timeBefore;
    trig_attributes{j}=Xa(find(~isnan(Xa(:,j))),j);
    C(j)=numel(trig_events{j});
end

[h bins]=hist(cell2mat(trig_events'),[timeBefore:binsize:timeAfter]);
