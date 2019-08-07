function [Iboot p5 p50 p95]= boostrap_Info_1D(f,sts,sts_prior,N,numBins)
%bootstrap repetitions of 1D information calculation for spike-triggered
%and prior distributions
%f- filter
%sts- spike triggered stimulus distribution
%stim - raw stimulus
%N number of bootstrap repetitions
%R number of draws from raw stimulus
%bins- number of bins for calculation
Iboot=nan(N,1);
n=size(sts,1);
parfor i=1:N
   Iboot(i)=MID_info(f,sts_prior,sts(randsample(n,n,true),:),numBins);
   if mod(i,5)==0,display(i/N);end
end


[H bins]=hist(Iboot,N*2);
H=H/sum(H);
c=cumsum(H);
p5=bins(find(c>.05));p5=p5(1);
p50=bins(find(c>.5));p50=p50(1);
p95=bins(find(c>.95));p95=p95(1);

