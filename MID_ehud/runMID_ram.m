function [t,MID,Ibest,I,sta]=runMID_ram(signal,fs,spike_times,timeBefore,timeAfter,staguess)
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
T_back=timeBefore; % seconds back (prior to spike)
T_for=timeAfter;  % seconds forward (after spike)

%check for proper dimensions %RAM
if size(spike_times,2)>size(spike_times,1), spike_times=spike_times';end
%sim. annealing definitions:
T0=30;       %initial temperature  %initial 1
dT=0.025;    %temperature step  %default .05
kmax=80;    %steps in each line maximization default 25
alpha=1;    %gradient step
Tmin=0.7;   %final temperature

%histogram binning
hist_bins=25;   %default 25
%% generate matrices
K_back=round(T_back*fs);
K_for=round(T_for*fs);
%spike indices
%round to make indices
ind_sp=round(spike_times*fs);
%kick out spikes that are close to edge of stim
ind_sp=ind_sp(ind_sp>K_back & ind_sp<=(length(signal)-K_for));

% %spike-triggered
% X=TriggeredSegments(signal,ind_sp,-K_for, K_back);

%spike triggering stimuli!
idx_sp=((((1-K_back):K_for)'*ones(1,length(ind_sp)))+(ones(K_back+K_for,1)*ind_sp'))'; %prior?
st_sp=signal(idx_sp)-mean(signal); %mean removed STA;
%sta guess for start
sta=mean(st_sp)';

%input indices
ind_stim=(K_back+1):(length(signal)-K_for);
%sampling prior (now hardcoded to 50000) samples RAM
%idx=((((1-K_back):K_for)'*ones(1,length(ind_stim)))+(ones(K_back+K_for,1)*ind_stim))'; %prior???

%corrected random samples??? ram
samps=randsample(numel(ind_stim),60000);
samps=samps(samps>K_back & samps<=(length(signal)-K_for));
samps=samps(1:50000);
idx=((((1-K_back):K_for)'*ones(1,50000))+(ones(K_back+K_for,1)*ind_stim(samps)))'; %prior???
st=signal(idx)-mean(signal);

% 
% prior
% X=TriggeredSegments(x,t,timeBefore, timeAfter);
figure
plot(sta)

%% run MID simulated annealing

if staguess==1
    [MID,Ibest,I]=MID_anneal_ram(st,st_sp,hist_bins,T0,dT,kmax,alpha,Tmin,sta);  %it calculates sta for start point
elseif numel(staguess)>1
    [MID,Ibest,I]=MID_anneal_ram(st,st_sp,hist_bins,T0,dT,kmax,alpha,Tmin,staguess);  %you can provide a guess
else
    [MID,Ibest,I]=MID_anneal_ram(st,st_sp,hist_bins,T0,dT,kmax,alpha,Tmin,0);   %no guess at all
end

if(mean(st_sp*MID)<0)
    MID=-MID;
end
t=[(1-K_back):K_for]/fs;






