function [GvI]=MID_grad(vector,stim,stim_sp,Nbins)
%   find Maximally Informative Dimentions (Sharpee et al)
%   Gradient ascent function
%   Code written by Avner Wallach
%   Do Not Distribute without authorization
%   % compute gradient of mutual information with respect to vector
% N- number of dimmentions in stimuli
% M- number of stims presented
% S- number of spikes detected
%
% inputs:
%       vector-      Nx1 vector describing direction of subspace
%       stim-           MxN matrix of stimuli samples
%       stim_sp-    SxN matrix of stimuli samples associated with spike
%       Nbins-           number of bins used in histogram calculations
% outputs:
%       GvI-               Nx1 gradient

% project to subspace
x=stim*vector; %Mx1
x_sp=stim_sp*vector;%Sx1

bins=linspace(min(x),max(x),Nbins+1);
Hv=histc(x,bins);
Pv=Hv(1:end-1)/sum(Hv);
Hv_sp=histc(x_sp,bins);
Pv_sp=Hv_sp(1:end-1)/sum(Hv_sp);
ddx=diff(Pv_sp./Pv);
GvI=0;
for i=1:Nbins-1
    ind=find(x>=bins(i+1) & x<bins(i+2));
    ind_sp=find(x_sp>=bins(i+1) & x_sp<bins(i+2));
    if(numel(ind_sp))
        GvI=GvI+(Pv(i+1)*(mean(stim_sp(ind_sp,:),1)-mean(stim(ind,:),1))'*ddx(i));
    end
end