function [ P, Tblock, c ,loc, cNdx] = PIB( T, K, maxiters, tau, beta )
%Takes in a transition probability matrix and outputs the permutation that 
%approximately block diagonalizes it using the predictive information 
%bottleneck algorithm, i.e. T(P,P) should be block diagonal
%
%   T(i,j) is pr transition from j to i.
%   K is the number of blocks
%   maxiters ~ maximum number of iterations
%   tau ~ number of time steps to consider
%   beta ~ complexity parameter (10 is good)

    if(~exist('beta','var'))
        beta=10;
    end
    if(~exist('tau','var'))
        tau=1;
    end

    N=size(T,1);

    przx=rand(K,N);
    przx=bsxfun(@times,przx,1./sum(przx,1));

    prx=mean((T^5000)')';
    prz=ones(K,1)/K;

    pryx=T^tau;
    pryz=pryx*bsxfun(@times,przx,prx')';


    for iters=1:maxiters

        DKL = sum(pryx.*log(pryx),1) - log(pryz')*pryx;
        lnprzx = bsxfun(@plus,-beta*DKL,log(prz));    
        lnprzx = bsxfun(@plus,lnprzx, -max(lnprzx));
        przx = bsxfun(@times,exp(lnprzx),1./sum(exp(lnprzx),1));
        prz=przx*prx;
        pryz=pryx*bsxfun(@times,przx,prx')';
        pryz=bsxfun(@times,pryz,1./sum(pryz,1));

    end


    [m,loc]=max(przx);
    c=unique(loc);
    
    for k=1:length(c)
        NA(k)=sum(loc==c(k));
    end
    [m,s]=sort(NA,'descend');
    c=c(s);  %unique labels sorted by size
    
    P=[];
    cNdx = cell(1,10);
    for k=1:length(c)
        P=[P,find(loc==c(k))];
        cNdx{k} = find(loc==c(k));
    end
    
    Tblock = T(P,P);
    
end

