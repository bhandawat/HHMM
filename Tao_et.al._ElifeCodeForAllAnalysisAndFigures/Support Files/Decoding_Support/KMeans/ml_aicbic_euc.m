function [aic,bic] = ml_aicbic_euc(data,cluster)
% [AIC, BIC] = ML_AICBIC_EUC(DATA,CLUSTERNUM)
% Calculate the Akaike Information Criterion and Baysian Information Criterion
% for a particular partitioning of multivariate data using an Euclidean
% distance measure
% data - a m*n feature matrix with m observations of n features
% clusternum - a k-element cell array.  Each containing the indices of
%              observations belonging to this cluster
% aic - returned AIC value
% bic - returned BIC value
%
% Xiang Chen, Dec 15, 2003

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

r = length(cluster);

for i=1:r

  idx=cluster{i} ;
  thiscluster = data(idx,:) ;

      clusterinfo(i).number = length(idx) ;
      clusterinfo(i).lambda = length(idx)/size(data, 1);%length(clusternum) ;
      if length(idx)==1
        clusterinfo(i).sumd = 0;
      else
        distances = ml_dist2(thiscluster,mean(thiscluster)) ;
        clusterinfo(i).sumd = sum(distances);
      end
end

lin = 0 ;
lin1 = 0 ;
lin2 = 0 ;
lin3 = 0 ;

for i=1:length(cluster)
    lin1 = lin1 + clusterinfo(i).number * log(clusterinfo(i).lambda) ;
    lin3 = lin3 + clusterinfo(i).sumd/2;
end

lin = lin1-lin2-lin3 ;
%
% the AIC as defined by Akaike is -2*lin + 2*k
%   where lin is the log of the maximum likelihood (calculated above), and
%   where k is the number of independently adjusted parameters
%
% for a Euclidean distance-based measure of likelihood,
%   (r-1) for the estimates of the membership of each cluster
%   n*r for the means for each variable for each cluster
%
Np = r-1+r*n;
aic = -2*lin + 2*Np;
bic = real(-2 * lin + Np * log(n));