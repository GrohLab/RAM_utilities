
function M = squishcell(C)

m=cellfun(@(x) size(x,1),C)

M=nan(max(m),numel(C));

for i=1:numel(C)
    M(1:numel(C{i}),i)=C{i};
end