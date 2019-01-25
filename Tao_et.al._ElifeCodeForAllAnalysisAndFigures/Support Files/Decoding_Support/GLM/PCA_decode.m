function [regressors,PCAStats,HLScomp,reconstructed] = PCA_decode(input)
    PCAStats = [];
    ndx = find(var(input)==0);
    input(:,ndx)=[];
    %w = 1./var(input);
    
    [coeff,score,latent,tsquared,explained,mu] = pca(input);
    PCAStats.wcoeff = coeff;
    PCAStats.score = score;
    PCAStats.latent = latent;
    PCAStats.tsquared = tsquared;
    PCAStats.explained = explained;
    
    % use this if eliminating PC based on percent explained
%     cdfExplained = zeros(1,length(explained));
%     for i = 1:length(explained)
%         cdfExplained(i+1) = cdfExplained(i)+explained(i);
%     end
%     idx = [1:min(find(cdfExplained>80))-1];

    % use this if eliminating PC based on eigenvalues
    percExplained = cumsum(latent)./sum(latent);
    idx = find(percExplained>0.9);
    idx = 1:min(idx);
    
    regressors = score(:,idx);
    HLScomp = [coeff(:,idx);latent(idx)'];
    c = HLScomp;
    for i = 1:length(ndx)
        c = insertrows(c, 0,ndx(i)-1);
    end
    HLScomp = c;
    
    reconstructed = score(:,idx) * coeff(:,idx)' + repmat(mu, size(input,1), 1);
    c = reconstructed;
    for i = 1:length(ndx)
    c = insertrows(c.', 0,ndx(i)-1).';
    end
    reconstructed = c;
end