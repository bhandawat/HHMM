function [aic,bic] = ml_aicbic_mah(data,clusternum, covs)
% [AIC, BIC] = ML_AICBIC_MAH(DATA,CLUSTERNUM, COVS)
% Calculate the Akaike Information Criterion and Baysian Information Criterion
% for a particular partitioning of multivariate data using a Mahalanobis
% distance measure on a global covariance matrix.
% data - a m*n feature matrix with m observations of n features
% clusternum - a k-element cell array.  Each containing the indices of
%              observations belonging to this cluster
% covs - an optional input of the covariance matrix used in distance
%        calculation.  It will be calculated from the whole data if omitted
% aic - returned AIC value
% bic - returned BIC value
%
% M. Boland - 29 Apr 1999
% M. VElliste, July 13, 2001. Removed debug printouts
% R.F. Murphy, July 17, 2001. Reinsert code to check bad cov
%                             Add check for not enough events in cluster
%                             Fix lin3 sum.  Change number of "free params"
%                             calculation in final aic calc to 2*n*r.
%                             Remove unnecessary "centers" argument.
%                             Add break on insufficient events in cluster.
% R.F. Murphy, July 19, 2001. Ignore small clusters in AIC calc.
%                             Return "free params" calc to Np
%                             Add documentation
% R.F. Murphy, July 21, 2001. Use rm_mahal to handle singular cov matrices.
% X. Chen, Aug 15, 2003. Using global covariance matrix.  Output BIC

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu


n=size(data,2) ;
r = length(clusternum);
if (~exist('covs', 'var'))
   covs = cov(data);
end

cluster = [];
for i=1:r
    cluster(clusternum{i}) = i; %[cluster clusternum{i}];
end
aic=0 ;

for i=1:r
    idx=find(cluster==i) ;
    thiscluster = data(idx,:) ;

    clusterinfo(i).cov = cov(thiscluster) ;
    clusterinfo(i).number = length(idx) ;
    clusterinfo(i).mean = mean(thiscluster, 1);
    clusterinfo(i).lambda = length(idx)/length(cluster) ;
    
     
    clusterinfo(i).sumd=0 ;
    for j=1:clusterinfo(i).number
      clusterinfo(i).sumd = clusterinfo(i).sumd + (thiscluster(j,:) - ...
             clusterinfo(i).mean) * inv(covs) * ...
             ((thiscluster(j,:)-clusterinfo(i).mean))' ;
    end
end

if (~isnan(aic))
    lin = 0 ;
    lin1 = 0 ;
    lin2 = 0 ;
    lin3 = 0 ;

    for i=1:max(cluster)
      lin1 = lin1 + clusterinfo(i).number * log(clusterinfo(i).lambda) ;
      lin3 = lin3 + clusterinfo(i).sumd/2 ;
    end

    lin = lin1-lin2-lin3 ;
    % For Mahalanobis distance measure using global Cov matrix
    aic = real(-2 * lin + 2 * ( r - 1 + r * n + n * (n + 1) * 0.5));
    bic = real(-2 * lin + (r - 1 + r * n + n * (n + 1) * 0.5) * log(n));
end