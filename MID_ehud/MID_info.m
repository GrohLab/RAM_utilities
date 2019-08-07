function [I]=MID_info(vector,stim,stim_sp,Nbins);
%   Code written by Avner Wallach
%   Do Not Distribute without authorization
%   % compute mutual information between spike arrival and projections of
% stimulus space onto vector
% N- number of dimmentions in stimuli
% M- number of stimi presented
% S- number of spikes detected
%
% inputs:
%       vector-      Nx1 vector describing direction of subspace
%       stim-           MxN matrix of stimuli samples
%       stim_sp-    SxN matrix of stimuli samples associated with spike
%       Nbins-           number of bins used in histogram calculations
% outputs:
%       I-                  scalar mutual information

% project to subspace
x=stim*vector; %Mx1
x_sp=stim_sp*vector;%Sx1

bins=linspace(min(x),max(x),Nbins+1);
Hv=histc(x,bins);
Pv=Hv(1:end-1)/sum(Hv);
Hv_sp=histc(x_sp,bins);
Pv_sp=Hv_sp(1:end-1)/sum(Hv_sp);
I=nansum(Pv_sp.*log2(Pv_sp./Pv));