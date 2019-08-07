function KLcorrected =KullbackLeibler2(p,q,bins)


% 
% p=p/mean(p);q=q/mean(q);
% 
% temp=[p(:);q(:)];
%[temp bins]=hist(temp,numBins);

%%
binsize=bins(2)-bins(1);

P=hist(p,bins);P=P/sum(P*binsize);
P(find(P==0))=eps;
Q=hist(q,bins);Q=Q/sum(Q*binsize);
Q(find(Q==0))=eps;
KL=(P.*log2(P./Q)*binsize);

KL=KL(find(isfinite(KL)));
KLraw=sum(KL);

KLcorrected=KLraw;
    