

%load Ross data from VPL
load 'Z:\Ross\Experiments\Counts_by_Condition\M168_Counts_by_Condition.mat'
ConditionsExport=ConditionsExport(1:4);

mins=[];
for i=1:ncond
    mins=[mins;cellfun(@size,ConditionsExport(i).Counts, 'uniformoutput',0)];
end
mins=cell2mat(mins);ntrial=min(min(mins(:,2:2:end)));
% we want to use this for multiple comparisons for each neuron (non
% parametric 1 factor ANOVA for each neuron
dim=size(cell2mat(ConditionsExport(2).Counts(:,2)));
nneuron=dim(1);
ncond=4; %number experimental conditions
%fix for unequal trial counts!


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



%%

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



%%
 %[p,tbl,stats] = kruskalwallis(counts)

 figure
 [h bins]=hist(deltas,linspace(min(deltas(:)),max(deltas(:)),30));plot(bins,h)
 title('difference between experimental and spont count medians, pooled neurons')

 
 %%
 fractions=nan(3,ncond);%no sig change, smaller, bigger
 for i=1:ncond
     i
     fractions(1,i)=sum(H(:,i)==0); %not sig
     fractions(2,i)=sum(deltas(:,i)<=0  & H(:,i)==1); % sig decrease
     fractions(3,i)=sum(deltas(:,i)>0 & H(:,i)==1); %sig increase
 end

 figure
 
 for i=1:ncond
    subplot(2,2,i)
    pie(fractions(:,i))
    title(ConditionsExport(i).name)
 end
 legend('not significant p>.05 Wilcoxon sign-rank','significant decrease','significant increase')
 %
 fractions=fractions/nneuron;
 
 
%% 
Spont=reshape(spont,size(spont,1),size(spont,2)*size(spont,3))
% counts1=cell2mat(cellfun(@(x) median(x'),ConditionsExport(2).Counts','UniformOutput',false))';
% counts2=cell2mat(cellfun(@(x) median(x'),ConditionsExport(3).Counts','UniformOutput',false))';

% counts=[counts1 counts2];

counts=[ squeeze(median(exp,2)) median(Spont')'];
%%volin plot practice
figure
violins = violinplot(counts, 'ShowNotches')
% spont=[violins(1).ScatterPlot.XData' violins(1).ScatterPlot.YData']
% exp=[violins(2).ScatterPlot.XData' violins(2).ScatterPlot.YData']
% 
% for i=1:size(spont,1)
%     plot([spont(i,1) exp(i,1)],[spont(i,2) exp(i,2)],'color',[.7 .7 .7],'linewidth',.1);
%     hold on
% end



%% normalize to relative change
counts=cell2mat(cellfun(@(x) median(x'),ConditionsExport(3).Counts','UniformOutput',false))';
delta=counts(:,2)./counts(:,1);
figure
violins = violinplot(delta(isfinite(delta)), 'ShowNotches')
median(delta(isfinite(delta)))
%%
P = signrank(counts(:,1),counts(:,2)) 
P = ranksum(counts(:,1),counts(:,3)) 
P = signrank(counts(:,2),counts(:,4)) 


 
 %% 
%  Table 2 Pairwise post-hoc Wilcoxon signed-rank tests of the mean EMG activites (Baseline (BL), 20%, 50%, 80% MVC) of the involuntarily contracting FDIpMA (z?=?z-value of Wilcoxon signed-rank test, p-value Dunn-Bonferroni adjusted for multiple comparisons, r?=?Pearson effect size, *significance).
% From: Structural Neural Correlates of Physiological Mirror Activity During Isometric Contractions of Non-Dominant Hand Muscles

%% Multiple comparisons (non-paired): control, L6 control, mech control, 1Hz L6, 10 Hz L6





