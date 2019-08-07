function C=colors(Cmap,num)

eval(['colormap ' Cmap])
cmap=colormap;

x=floor(63/num);
indices=[1:x:(x*num)];
C=cmap(indices,:);