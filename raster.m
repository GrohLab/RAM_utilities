

function r=raster(spikes,level,col,h)

%this function produces a raster plot of vertical lines for multiple
%spike trains, spikes, a structure;
%Each spike train is plotted at increasing vertical coordinates, with height h,
%in color specified by string col, e.g. 'k','r'
%r is the final handle produced, this can be used to turn the baseline off
%(deals with bizarre behavior of Matlab 2010 stem plot).

r=[];
if ~isempty(spikes)
    ys=	(ones(1,numel(spikes))*(level+h*.9));
    r=stem(spikes,ys,'basevalue',level);
    set(r,'color', col, 'marker','none','linewidth',1);
    set(r,'ShowBaseLine','on');
end


