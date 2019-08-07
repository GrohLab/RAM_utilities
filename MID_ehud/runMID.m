function [t,MID,Ibest,I]=runMID(signal,fs,spike_times)
%   find Maximally Informative Dimentions (Sharpee et al)
%   Code written by Avner Wallach
%   Do Not Distribute without authorization
%   
%inputs:
%   signal- continuous input signal (vector)
%   fs-     sampling frequency of signal (Hz)
%   spike_times- time of spike detections (in s)
%
%outputs:
%   t-   time vector for plotting the MID projection
%   MID- the MID projection
%   Ibest- mutual information in the MID projection
%   I- vector of mutual information calculated in MID iterations.
%
%   to see the convergence of the algorithm, plot(I)
%   to see the MID projection, plot(t,MID)
%
%% parameters
%timeframe definitions: 
T_back=0.5; % seconds back (prior to spike)
T_for=0.1;  % seconds forward (after spike)

%sim. annealing definitions:
T0=1;       %initial temperature
dT=0.05;    %temperature step
kmax=25;    %steps in each line maximization
alpha=1;    %gradient step
Tmin=0.7;   %final temperature

%histogram binning 
hist_bins=25;
%% generate matrices
K_back=round(T_back*fs);
K_for=round(T_for*fs);
%spike indices
ind_sp=round(spike_times*fs);
ind_sp=ind_sp(ind_sp>K_back & ind_sp<=(length(signal)-K_for));    
idx_sp=((((1-K_back):K_for)'*ones(1,length(ind_sp)))+(ones(K_back+K_for,1)*ind_sp'))';
st_sp=signal(idx_sp)-mean(signal);

%input indices
ind_stim=(K_back+1):(length(signal)-K_for);
idx=((((1-K_back):K_for)'*ones(1,length(ind_stim)))+(ones(K_back+K_for,1)*ind_stim))';
st=signal(idx)-mean(signal);

%% run MID simulated annealing
[MID,Ibest,I]=MID_anneal(st,st_sp,hist_bins,T0,dT,kmax,alpha,Tmin);
if(mean(st_sp*MID)<0)
    MID=-MID;
end
t=[(1-K_back):K_for]/fs;