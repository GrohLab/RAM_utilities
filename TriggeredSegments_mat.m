function X=TriggeredSegments_mat(M,t,timeBefore, timeAfter)
timeBefore=round(timeBefore);timeAfter=round(timeAfter);

Xmat=spalloc(-timeBefore+timeAfter+1,numel(t),);
for i=1:numel(t)
    indices=(t(i)+timeBefore):(t(i)+timeAfter);
    % Omitting the segment of the signal if the resquested time exceeds the
    % actual domain of the signal.
    if indices(1) < 1
        X(:,i) = NaN;
    else
        X(:,i)=x(round(indices));
    end
end
