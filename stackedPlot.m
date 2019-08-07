function f=stackedPlot(x,Y,order,spacer)

if size(Y,2)~=numel(x),Y=Y';end
if size(Y,2)~=size(x,2),x=x';end

if order ~=0
Y=Y(order,:);
end

m=[0 max(Y)];
f=figure

for i=1:size(Y,1)
    plot(x,Y(i,:)+m(i)*spacer*i,'k')
    hold on
end


