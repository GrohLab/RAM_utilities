function [h bins trig_events C]=eventHistogram(events,trigs,N,binsize,timeBefore,timeAfter)


bievents=zeros(N,1);
bievents(events)=1;


X=TriggeredSegments(bievents,trigs,timeBefore,timeAfter);
trig_events={};
C=[];
for j=1:(size(X,2))
    trig_events{j}=find(X(:,j))+timeBefore;
    C(j)=numel(trig_events{j});
end
C=C';
[h bins]=hist(cell2mat(trig_events'),[timeBefore:binsize:timeAfter]);

