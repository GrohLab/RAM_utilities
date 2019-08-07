cd C:\Users\rebecca\Dropbox\SynchronyFall2014
clear all
%load all 23 single cells
%script
GetSinglePOm_cells
DS=20; %downsampling for eeg preprocessing (to 1000 Hz)
%%
MID={};
Ibest={};
Inf={};
stas={}
%%
cd C:\Users\rebecca\Dropbox\synchronyFall2014\MID_ehud
timeBefore=.4
timeAfter=0;
for I=1:23
    ds=1000;
    spikes=Data{I}.spikes/20000;
    eeg=Data{I}.eeg_clean; eeg=zscore(eeg);
    fs=500
    tic
    I
    [t,MID{I},Ibest{I},Inf{I},stas{I}]=runMID_ram(eeg(2:2:end),fs,spikes(1:1:end),timeBefore,timeAfter,0);
    [t,MID_sta{I},Ibest_sta{I},Inf_sta{I},stas{I}]=runMID_ram(eeg(2:2:end),fs,spikes(1:1:end),timeBefore,timeAfter,1);
    toc
end

save POM_mid_spikes_t10 MID Ibest Inf stas MID_sta Ibest_sta Inf_sta


%%  run events only with new parameters, parallelized for speed
pack
MID={};
Ibest={};
Inf={};
stas={}
FirstSpikes={};
cd C:\Users\rebecca\Dropbox\synchronyFall2014\MID_ehud
timeBefore=.4
timeAfter=0;
ds=1000;
%grab variables for parallel run
fs=500
for I=1:23
    spikes=Data{I}.spikes;
    %extract bursts
    [Bursts firstSpikes Cs{I}]=returnBursts(spikes,20*20);
    spikes=firstSpikes/20000; %first spikes of events, in seconds
    eeg=Data{I}.eeg_clean; eeg=zscore(eeg);
    EEGs{I}=eeg(2:2:end); %downsample to 500 Hz
    FirstSpikes{I}=spikes;
end



%%

parfor I=1:23
    display(I)
    [t,MID{I},Ibest{I},Inf{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,0);
    [t,MID_sta{I},Ibest_sta{I},Inf_sta{I},stas{I}]=runMID_ram(EEGs{I},fs,FirstSpikes{I},timeBefore,timeAfter,1);
    %     save POM_mid_events_t6 MID Ibest Inf stas MID_sta Ibest_sta Inf_sta Cs FirstSpikes t
end
save POM_mid_events_t6 MID Ibest Inf stas MID_sta Ibest_sta Inf_sta Cs FirstSpikes



%%
plot(t,zscore(MID{1}),t,zscore(MID_sta{1}))
figure
for I=1:numel(Data)
    mid{I}=MID{I}/std(MID{I});%-mean(MID{I});
    plot(t,mid{I})
    pause
end
plot(cell2mat(mid))


%% compare features
Data={};
events=load('POM_mid_events')

Data{1}=load('POM_mid_spikes')
Data{2}=load('POM_mid_spikes_t10')


Data{1}=load('POM_mid_events')
Data{2}=load('POM_mid_events_t6')
%%
figure
for i=1:23
    t=[-200:-1]/fs;
    subplot(1,2,1)
    f1=Data{1}.MID{i};
    f2=Data{2}.MID{i};
    f3=Data{1}.MID_sta{i};
    f4=Data{2}.MID_sta{i};
    %f4=Data{1}.stas{i};
    f1=f1/norm(f1);
    f2=f2/norm(f2);
    f3=f3/norm(f3);
    f4=f4/norm(f4)
    p=plot(t,f1,t,f2,t,f3,t,f4)
    set(p,'linewidth',2)
    legend('t=1','t=6','sta seeded','STA')
    
    hold off
    subplot(1,2,2)
    I1=Data{1}.Inf{i};
    I2=Data{2}.Inf{i};
    I3=Data{1}.Inf_sta{i};
    I4=Data{2}.Inf_sta{i};
    
    plot(I1)
    hold on
    plot(I2);
    plot(I3);
    plot(I4);
    legend('t=1','t=10','sta seeded t=1', 'sta seeded t=6')
    pause
    hold off
end



%% check different event sizes
load 'POM_mid_events_t6' % with t=6 starting temp
 
MID0=MID;
load newMID_events
MID1=MID
load newMID_events2
MID2=MID;

figure
for I=1:numel(MID)
    plot(MID0{I})
    hold on
       plot(MID1{I})
        plot(MID2{I})
        title(I)
        pause
    hold off
end
%%  REDO:  16 22
load 'POM_mid_events_t6' % with t=6 starting temp
load BestEventMID_sofar MID

close all
Info_mid=[];
Info_sta=[];
Info1=[];
Info2=[];
BP=[];
for I=[20:23]%sort([22 21 18 17 16 12 11 10])  %1:23
    I
    fs=500;
    t=[-200:-1]/fs;
    f=MID_sta{I};f=f/norm(f);
    stim=EEGs{I};
    events=FirstSpikes{I}*20000; %first spikes in indices for 20000 kHz
    events=round(events/(20000/fs));
    goods=find(events>numel(t) & events<numel(stim));  %kick out bad events
    events=events(goods);C=Cs{I}(goods); %adjust counts vector to maintain indexing
    sts=TriggeredSegments(stim,events,-(numel(t)-1), 0)'; %get spike-triggered stim;
    STA=mean(sts)';STA=STA/norm(STA);
    prior_samp=randsample(numel(stim),201000);
    prior_samp=prior_samp(find(prior_samp>numel(t)));prior_samp=prior_samp(1:200000);
    sts_prior=TriggeredSegments(stim,prior_samp,-(numel(t)-1), 0)';
    %sts_prior=sts;
    Info_mid(I)=MID_info(f,sts_prior,sts,50);
    Info_sta(I)=MID_info(STA,sts_prior,sts,50);
    figure
    subplot(2,2,1);plot(t,STA);subplot(2,2,2),plot(t,f)
    % now check by event size
    subplot(2,2,3)
    
    hist(C);
    title(I)
    c1=find(C==1);     % 12 16 22 21
    c2=find(C>1);
    
    sts1=TriggeredSegments(stim,events(c1),-(numel(t)-1), 0)'; %get spike-triggered stim;
    sts2=TriggeredSegments(stim,events(c2),-(numel(t)-1), 0)'; %get spike-triggered stim;
    Info1(I)=MID_info(f,sts_prior,sts1,50);
    if numel(c2)>10
        Info2(I)=MID_info(f,sts_prior,sts2,50);
    else
        Info2(I)=nan;
    end
    BP(I)=numel(c2/numel(events));
    
    subplot(2,2,4)
    bar([1:4],[Info_sta(I) Info_mid(I) Info1(I) Info2(I)])
    pause
    
end


%% 17 9 6 5 4 3 1  Bursty examples
figure
plot(Info1,Info2,'o')
hold on
plot([0:4],[0:4])


%%
counter=0
close all
figure

info1=[]
info2=[]
info3=[]

for I=[17 9 6 5 4 3 1 ]
    
    
    I
    fs=500;
    t=[-200:-1]/fs;
    f=MID{I};f=f/norm(f);
    stim=EEGs{I};
    events=FirstSpikes{I}*20000; %first spikes in indices for 20000 kHz
    events=round(events/(20000/fs));
    goods=find(events>numel(t) & events<numel(stim));  %kick out bad events
    events=events(goods);C=Cs{I}(goods); %adjust counts vector to maintain indexing
    sts=TriggeredSegments(stim,events,-(numel(t)-1), 0)'; %get spike-triggered stim;
    STA=mean(sts)';STA=STA/norm(STA);
    prior_samp=randsample(numel(stim),101000);
    prior_samp=prior_samp(find(prior_samp>numel(t)));prior_samp=prior_samp(1:100000);
    sts_prior=TriggeredSegments(stim,prior_samp,-(numel(t)-1), 0)';
    %sts_prior=sts;
    
    c1=find(C==1);
    c2=find(C==2);
    c3=find(C==3);
    c4=find(C>3);
    
    
    sts1=TriggeredSegments(stim,events(c1),-(numel(t)-1), 0)'; %get spike-triggered stim;
    sts2=TriggeredSegments(stim,events(c2),-(numel(t)-1), 0)'; %get spike-triggered stim;
    sts3=TriggeredSegments(stim,events(c3),-(numel(t)-1), 0)'; %get spike-triggered stim;
    sts4=TriggeredSegments(stim,events(c4),-(numel(t)-1), 0)'; %get spike-triggered stim;
    
    counter=counter+1
    
    info1(counter)=MID_info(f,sts_prior,sts1,50);
    info2(counter)=MID_info(f,sts_prior,sts2,50);
    info3(counter)=MID_info(f,sts_prior,sts3,50);
    info4(counter)=MID_info(f,sts_prior,sts4,50);
    
    subplot(2,4,counter)
    bar([1:4],[info1(counter) info2(counter)  info3(counter) info4(counter)])
    
    title(I)
end


%%

figure
bar(mean([info1' info2' info3' info4']))
bar(mean([info1'./info1' info2'./info1' info3'./info1']))


%% bootstrapping error, STA and MID feature, for synchrony paper

close all
Info_mid=[];
Info_sta=[];

Iboot_mid={}
p5=[];p50=[];p95=[];
Iboot_sta={}
p5_sta=[];p50_sta=[];p95_sta=[];
GOODS=[1:11 13:15 17:23]
for I=GOODS
    I
    fs=500;
    t=[-200:-1]/fs;
    f=MID{I};f=f/norm(f);
    stim=EEGs{I};
    events=FirstSpikes{I}*20000; %first spikes in indices for 20000 kHz
    events=round(events/(20000/fs));
    goods=find(events>numel(t) & events<numel(stim));  %kick out bad events
    events=events(goods);C=Cs{I}(goods); %adjust counts vector to maintain indexing
    sts=TriggeredSegments(stim,events,-(numel(t)-1), 0)'; %get spike-triggered stim;
    STA=mean(sts)';STA=STA/norm(STA);
    prior_samp=randsample(numel(stim),201000);
    prior_samp=prior_samp(find(prior_samp>numel(t)));prior_samp=prior_samp(1:200000);
    sts_prior=TriggeredSegments(stim,prior_samp,-(numel(t)-1), 0)';
    N=500;
    [Iboot_mid{I} p5(I) p50(I) p95(I)]= boostrap_Info_1D(f,sts,sts_prior,N,50);
    [Iboot_sta{I} p5_sta(I) p50_sta(I) p95_sta(I)]= boostrap_Info_1D(STA,sts,sts_prior,N,50);
    Info_mid(I)=MID_info(f,sts_prior,sts,50);
   Info_sta(I)=MID_info(STA,sts_prior,sts,50);
end


%% calculate  p values, plot stas and mid features
mid=[]
sta=[];
for i=GOODS
    temp=stas{i};
    temp=zscore(temp);temp=temp/norm(temp);
    sta=[sta; temp'];
    
    
    temp=MID{i};
    temp=zscore(temp);
    temp=temp/norm(temp)
    mid=[mid; temp'];
end

m=mean(sta);s=std(sta);
figure
plot(t,m,t,m-s,t,m+s)
hold on
m=mean(mid);s=std(mid);
plot(t,m,t,m-s,t,m+s)



%% plot examples of STA and MID and projections on to each
H_sta_sp=[];
H_sta_prior=[];
H_mid_sp=[];
H_mid_prior=[];
for I=GOODS
    sta=stas{I};
    sta=zscore(sta);sta=sta/norm(sta);
    mid=MID{I};
    mid=zscore(mid);
    mid=mid/norm(mid)
    
    figure
    subplot(2,2,[1 3])
    plot(t,sta,t,mid)
    
 
    f=MID{I};f=f/norm(f);
    stim=EEGs{I};
    events=FirstSpikes{I}*20000; %first spikes in indices for 20000 kHz
    events=round(events/(20000/fs));
    goods=find(events>numel(t) & events<numel(stim));  %kick out bad events
    events=events(goods);C=Cs{I}(goods); %adjust counts vector to maintain indexing
    sts=TriggeredSegments(stim,events,-(numel(t)-1), 0)'; %get spike-triggered stim;
    
    STA=sta;STA=STA/norm(STA);
    prior_samp=randsample(numel(stim),201000);
    prior_samp=prior_samp(find(prior_samp>numel(t)));prior_samp=prior_samp(1:200000);
    sts_prior=TriggeredSegments(stim,prior_samp,-(numel(t)-1), 0)';
    
    s_sp_sta=STA'*sts';
    s_sp_mid=f'*sts';
    s_sta=STA'*sts_prior';
    s_mid=f'*sts_prior';
    
    s_sp_sta=s_sp_sta/std(s_sta)-mean(s_sta);
    s_sta=s_sta/std(s_sta)-mean(s_sta);

    s_sp_mid=s_sp_mid/std(s_mid)-mean(s_mid);
    s_mid=s_mid/std(s_mid)-mean(s_mid);
    
     [h bins]=hist([s_sp_sta s_sp_mid s_sta s_mid],[-10:.025:10]);
     H=histnorm({s_sp_sta' s_sta' s_sp_mid' s_mid'},bins);
     subplot(2,2,2)
     plot(bins,H(:,1:2))
     legend('spike-triggered sta','prior sta')
     subplot(2,2,4)
     plot(bins,H(:,3:4))
     legend('spike-triggered MID', 'prior mid')
    H_sta_sp=[H_sta_sp H(:,1)];
    H_sta_prior=[H_sta_prior H(:,2) ];
    H_mid_sp=[H_mid_sp H(:,3)];
    H_mid_prior=[H_mid_prior H(:,4)];
end


%%
close all
xlims=[-4 5]
ylims=[0 .05]
figure
subplot(2,1,1)

axis([xlims ylims])
m1=mean(H_sta_sp');
s1=std(H_sta_sp');
m=mean(H_sta_prior');
s=std(H_sta_prior');
plot(bins,m,bins,m1)%m-s,bins, m+s,bins,m1,bins,m1-s1,bins,m1+s1)
axis([xlims ylims]) 
subplot(2,1,2)
m1=mean(H_mid_sp');
s1=std(H_mid_sp');
m=mean(H_mid_prior');
s=std(H_mid_prior');
plot(bins,cumsum(m),bins,cumsum(m1%,m-s,bins, m+s,bins,m1,bins,m1-s1,bins,m1+s1)

axis([xlims ylims])
%% ploting information distributions from bootstrapping
figure
pval_info=[]
for i=GOODS
[h,pval_info(i)] = kstest2(Iboot_mid{i},Iboot_sta{i},'tail','smaller');   
% pval_info(i) = ranksum(Iboot_mid{i}',Iboot_sta{i}','alpha',.05)
hist([Iboot_mid{i},Iboot_sta{i}],50)
title(i)
pause
end


%% plot MID info vs STA info
figure
ploterr(p50_sta,p50,{p5_sta,p95_sta},{p5,p95},'k.')
hold on
plot([0 4.5],[0 4.5],':k')
axis([0 4.5 0 4.5 ])
axis square
xlabel ('Info_{ETA} [bits]')
ylabel ('Info_{MID} [bits]')


%% compare eeg to gaussian
figure
eeg_all=cell2mat(EEGs');s=std(eeg_all);m=mean(eeg_all);
[Heeg bins_eeg]=hist(eeg_all,1000);
[H bins]=histnorm(EEGs,bins_eeg);
%%
figure
G= exp((-(bins_eeg-m).^2)/(2*v))./((2*v*pi)^.5);G=G/sum(G);
plot(bins_eeg,G,'k:','linewidth',4)
hold on
plot(bins_eeg,H,'linewidth',1)
plot(bins_eeg,G,'k:','linewidth',4)
xlim([-6 6])
legend('Unit gaussian')
legend boxoff
xlabel 'Normalized LFP amplitude (z-score)'
ylabel P(LFP)
title 'Distribution of LFP is non-gaussian (n=23)'
p=[];
%stat test for normal distribution 
for i=1:numel(EEGs)
    [h,p(i)] = kstest(EEGs{i});
end

%%
figure
plot(Info_sta,Info_mid,'o')
hold on
plot([0 4],[0 4],'k:')
axis square
xlabel('STA information [bits]')
ylabel('MID information [bits]')




%% now information for burst and 