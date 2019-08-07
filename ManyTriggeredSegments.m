function Xs=ManyTriggeredSegments(X,t,timeBefore, timeAfter);
timeBefore=ceil(timeBefore);timeAfter=ceil(timeAfter);


if iscell(X)
    Xs={};
    for i=1:numel(X)
        Xs{i}=TriggeredSegments(X{i},t,-timeBefore, timeAfter);
    end
else
    Xs=TriggeredSegments(X{i},t,-timeBefore, timeAfter);
end