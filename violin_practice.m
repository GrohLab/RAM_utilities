clear all


ncond=4; %number experimental conditions
%fix for unequal trial counts!

%load Ross data from VPL


files={'Z:\Ross\Experiments\Counts_by_Condition\M168_Counts_by_Condition.mat'
    'Z:\Ross\Experiments\Counts_by_Condition\M160_Counts_by_Condition.mat'
    'Z:\Ross\Experiments\Counts_by_Condition\M167HL_Cortex_Counts_by_Condition.mat'
    'Z:\Ross\Experiments\Counts_by_Condition\M167VPL_Counts_by_Condition.mat'
    'Z:\Ross\Experiments\Counts_by_Condition\M168_Counts_by_Condition.mat'}
I=3;


[FILEPATH,NAME,EXT] = fileparts(files{I});
load(files{I})
ConditionsExport=ConditionsExport(1:4);


%%
dim=size(cell2mat(ConditionsExport(2).Counts(:,2)));
nneuron=dim(1);
ncond=4;

mins=[];
for i=1:ncond
    mins=[mins;cellfun(@size,ConditionsExport(i).Counts, 'uniformoutput',0)];
end
mins=cell2mat(mins);ntrial=min(min(mins(:,2:2:end)));
% we want to use this for multiple comparisons for each neuron (non
% parametric 1 factor ANOVA for each neuron

exp=nan(nneuron,ntrial,ncond);%preallocate;
spont=nan(nneuron,ntrial,ncond);
ccount=0;
time=1; %placeholder for time chunk to convert to rates;
for i=1:ncond
    ccount=ccount+1;
    newexp=cell2mat(ConditionsExport(i).Counts(:,2))/time;
    newspont=cell2mat(ConditionsExport(i).Counts(:,1))/time; %pool spontaneous.
    exp(:,:,ccount)=newexp(:,1:ntrial);
    spont(:,:,ccount)=newspont(:,1:ntrial);
end


%
%I should add something for probability-- we can do chi2 or binomial test.
%What about paired aspect?

%exp_bi
%spont_bi

% first rank-sum (paired by pairs of observations) for each experimental
% condition

p=nan(nneuron,ncond);
deltas=nan(nneuron,ncond);
for i=1:ncond
    for n=1:nneuron
        p(n,i) =  signrank(exp(n,:,i)',spont(n,:,i)');
        deltas(n,i)=median(exp(n,:,i)')-median(spont(n,:,i)');
    end
end
H=p;H(H>.05)=0;H(H~=0)=1;  %each column a different experimental condition

decrease=(deltas<0 & H==1)
increase=(deltas>0 & H==1)
nochange=H==0;

ClusterInfo.Spont.increase=increase;
ClusterInfo.Spont.decrease=decrease;
ClusterInfo.Spont.nochange=nochange;

% H is nneurons by ncond, tells us if experimental changes
%%[p,tbl,stats] = kruskalwallis(counts)

figure
[h bins]=hist(deltas,linspace(min(deltas(:)),max(deltas(:)),30));plot(bins,h)
title('difference between experimental and spont count medians, pooled neurons')


%
fractions=nan(3,ncond);%no sig change, smaller, bigger
for i=1:ncond
    i
    fractions(1,i)=numel(find(H(:,i)==0)); %not sig
    fractions(2,i)=numel(find(deltas(:,i)<0  & H(:,i)==1)); % sig decrease
    fractions(3,i)=sum(deltas(:,i)>0 & H(:,i)==1); %sig increase
end

figure



for i=1:ncond
    subplot(2,2,i)
    mypie=pie(fractions(:,i))
    
    %make color nicer!
    x=title([ConditionsExport(i).name ' ' NAME])
    set(x, 'Interpreter', 'none')
    
end

legend('not sig. p>.05 W sr','decrease','increase')


%  compare mech control to conditions with laser=========================



p=nan(nneuron,2);
deltas=nan(nneuron,2);
for i=1:2
    for n=1:nneuron
        p_laser(n,i) = ranksum(exp(n,:,i)',exp(n,:,4)');
        deltas_laser(n,i)=median(exp(n,:,i)')-median(exp(n,:,4)');
    end
end
Hlaser=p_laser;Hlaser(Hlaser>.05)=0;Hlaser(Hlaser~=0)=1;  %each column a different experimental condition

decrease=(deltas_laser<0 & Hlaser==1)
increase=(deltas_laser>0 & Hlaser==1)
nochange=Hlaser==0;


ClusterInfo.MechLight.increase=increase;
ClusterInfo.MechLight.decrease=decrease;
ClusterInfo.MechLight.nochange=nochange;


% cell breakdown
fractions=nan(3,2);%no sig change, smaller, bigger
for i=1:2
    i
    fractions(1,i)=numel(find(Hlaser(:,i)==0)); %not sig
    fractions(2,i)=numel(find(deltas_laser(:,i)<0  & Hlaser(:,i)==1)); % sig decrease
    fractions(3,i)=sum(deltas_laser(:,i)>0 & Hlaser(:,i)==1); %sig increase
end


figure
for I=1:2
    counts=[ squeeze(median(exp(:,:,4),2)) squeeze(median(exp(:,:,I),2))];
    %%volin plot practice
    subplot(2,3,[2 3]+(I-1)*3)
    violins = violinplot(counts, {'Mech', 'Mech+Laser'})
    spontpoints=[violins(1).ScatterPlot.XData' violins(1).ScatterPlot.YData']
    exppoints=[violins(2).ScatterPlot.XData' violins(2).ScatterPlot.YData']
    
    for i=1:size(exp,1)
        if Hlaser(i,I)
            if deltas_laser(i,I)>0
                plot([spontpoints(i,1) exppoints(i,1)],[spontpoints(i,2) exppoints(i,2)],'-','color',[.7 .7 .7]*0,'linewidth',1);
            elseif deltas_laser(i,I)<0
                plot([spontpoints(i,1) exppoints(i,1)],[spontpoints(i,2) exppoints(i,2)],'-','color',[.7 0 0],'linewidth',1);
            end
        end
        hold on
    end
    title(ConditionsExport(I).name)
    ylabel 'median counts'
end





for i=1:2
    subplot(2,3,(3*i-2))
    mypie=pie(fractions(:,i))
    
    %make color nicer!
    x=title([ConditionsExport(i).name ' ' NAME])
    set(x, 'Interpreter', 'none')
    
end

legend('not sig. p>.05 W sr','decrease','increase')


%% normalize to relative change
counts=cell2mat(cellfun(@(x) median(x'),ConditionsExport(3).Counts','UniformOutput',false))';
delta=counts(:,2)./counts(:,1);
figure
violins = violinplot(delta(isfinite(delta)), 'ShowNotches')
median(delta(isfinite(delta)))
%%



%%
%  Table 2 Pairwise post-hoc Wilcoxon signed-rank tests of the mean EMG activites (Baseline (BL), 20%, 50%, 80% MVC) of the involuntarily contracting FDIpMA (z?=?z-value of Wilcoxon signed-rank test, p-value Dunn-Bonferroni adjusted for multiple comparisons, r?=?Pearson effect size, *significance).
% From: Structural Neural Correlates of Physiological Mirror Activity During Isometric Contractions of Non-Dominant Hand Muscles

%% Multiple comparisons (non-paired): control, L6 control, mech control, 1Hz L6, 10 Hz L6




%%
%fractions=fractions/nneuron;


%
Spont=reshape(spont,size(spont,1),size(spont,2)*size(spont,3))
% counts1=cell2mat(cellfun(@(x) median(x'),ConditionsExport(2).Counts','UniformOutput',false))';
% counts2=cell2mat(cellfun(@(x) median(x'),ConditionsExport(3).Counts','UniformOutput',false))';
% counts=[counts1 counts2];
counts=[ squeeze(median(exp,2)) median(Spont')'];
%%volin plot practice
figure
violins = violinplot(counts, 'ShowNotches')
spont=[violins(1).ScatterPlot.XData' violins(1).ScatterPlot.YData']
exp=[violins(2).ScatterPlot.XData' violins(2).ScatterPlot.YData']

for i=1:size(spont,1)
    plot([spont(i,1) exp(i,1)],[spont(i,2) exp(i,2)],'color',[.7 .7 .7],'linewidth',.1);
    hold on
end





