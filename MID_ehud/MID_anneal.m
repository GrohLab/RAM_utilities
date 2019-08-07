function [vbest,Ibest,I]=MID_anneal(stim,stim_sp,Nbins,T0,dT,kmax,alpha,Tmin)
%   find Maximally Informative Dimentions (Sharpee et al)
%   simulated annealing function
%   Code written by Avner Wallach
%   Do Not Distribute without authorization
%
N=size(stim,2);
v(:,1)=rand(N,1); %initial guess
I(1)=MID_info(v(:,1),stim,stim_sp,Nbins); %initial information
vbest=v(:,1);
Ibest=I(1);
m=1;
T=T0;
while(T>Tmin)% line maximizations
    k=1;
    while(k<kmax)
        [GvI]=MID_grad(v(:,m),stim,stim_sp,Nbins); %compute gradient
        v_new=v(:,m)+alpha*GvI; %gradient ascent;
        I_new=MID_info(v_new,stim,stim_sp,Nbins);
        if(I_new>I(m) | exp((I_new-I(m))/T)>rand())
            m=m+1;
            v(:,m)=v_new;
            I(m)=I_new;
            if(I_new>Ibest)
                Ibest=I_new;
                vbest=v_new;
            end
        end
        k=k+1;
    end
    v(:,m)=vbest;   %return to best vector
    I(m)=Ibest;   %return to best info
    T=T*(1-dT); %decrease temperature
    alpha=alpha*(1-dT);
end
end
