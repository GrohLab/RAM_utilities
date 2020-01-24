function F=returnFields(S,Str)
if ~iscell(Str)
    F=returnFields_helper(S,Str)
else
% f={};
%     for j=1:numel(Str)
%         str=Str{j}
%         f{j}=returnFields_helper(S,str);
%     end
%     F=struct('temp',Str)
%     for j=1:numel(Str)
%        F.(Str{j})=f{j};
%     end
end




function f = returnFields_helper(S,Str)
f=cell(size(S));
for i=1:numel(S)
    f{i}=getfield(S{i},Str);
end