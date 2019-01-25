function [kMeans] = ml_kmeans_aicbic(features, resultfile, distance, idx, kmin, kmax, featidx, ...
    weight )
% FUNCTION ML_KMEANS_AICBIC(FEATURES, RESULTFILE, DISTANCE, FEATIDX, 
%   WEIGHT, KMIN, KMAX)
% Perform k-means algorithm with different k (from 2 to the number of clones in
% feat) followed by AIC and BIC calculation.  The optimal k could be get by
% minimizing AIC or BIC values.
% features: the input feature values. It is a cell array where each element is
%           the feature matrix for one clone
% resultfile: filename to store the result.  Use ml_readaicbic to read the
%             information
% distance: specify the distance function.  One of the following values: 
%           'euclidean' (by default), 'mahalanobis', 'both'.
% featidx: the index of features used in distance calculation.  If not
%          provided, all features will be included.
% weight: an optional input for weighting each feature.  Not applicable when
%         the method is 'mahalanobis'.  The weight is applied after all 
%         features are normalized
% Xiang Chen, Dec 22, 2004

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

if (~exist('distance', 'var'))
    distance = 'eu';
end

switch (distance(1:2))
case 'eu'
    eu = 1;
    ma = 0;
case 'ma'
    ma = 1;
    eu = 0;
case 'bo'
    eu = 1;
    ma = 1;
otherwise
    error('Distance must be one of the following ''euclidean'', ''mahlanobis'' or ''both''.  Type ''help ml_kmeans_aicbic'' for details.');
end

feat = [];
groupinfo = [];
% if (length(features) <= 20)
%      upper = 2;
% else
     upper = 1;
% end

for m = 1 : length(features)
     %load features into an array
     feat = [feat; double(features{m})];
     groupinfo = [groupinfo, m + (1:size(features{m}, 1)) / 1000];
end

if exist('featidx', 'var')
     feat = feat(:, featidx);
end

feat = zscore(feat);
if exist('weight', 'var')
     feat = feat .* repmat(weight, [size(feat, 1) 1]);
end

if (ma)
    cov_f = (cov(feat));
end

if ~exist('kmax','var')
    kmax = upper * length(features);
end

if ~exist('kmin','var')
    kmin = 2;
end

for m = kmin : kmax
    if (ma)
        [cluster_mahs{m}, group_mah{m-kmin+1}] = ml_kmeans_mah(feat, groupinfo, m,cov_f);
        [a, b] = ml_aicbic_mah(feat, cluster_mahs{m}, cov_f);
        aic_mah(m-kmin+1) = a;
        bic_mah(m-kmin+1) = b;
    end
    if (eu)
        [cluster_eucs{m}, group_euclid{m-kmin+1}] = ml_kmeans_euc(feat, groupinfo, m, idx);
        [a, b] = ml_aicbic_euc(feat, cluster_eucs{m});
        aic_euclid(m-kmin+1) = a;
        bic_euclid(m-kmin+1) = b;
    end
    %save('/imaging/tmp/tmp_cluster101_norange.mat', 'group_mah', 'group_euclid', 'aic_mah', 'bic_mah', 'aic_euclid', 'bic_euclid', 'cov_f', 'as_mah', 'as_euclid');
end

kMeans = [];
switch (distance(1:2))
case 'eu'
    kMeans.group = group_euclid;
    kMeans.cluster = cluster_eucs;
    kMeans.aic = aic_euclid;
    kMeans.bic = bic_euclid;
    kMeans.type = 'Euclidean';
case 'ma'
    kMeans.group = group_mah;
    kMeans.cluster = cluster_mahs;
    kMeans.aic = aic_mah;
    kMeans.bic = bic_mah;
    kMeans.cov_f = cov_f;
    kMeans.type = 'Mahalanobis';
case 'bo'
    kMeans.euc.group = group_euclid;
    kMeans.euc.cluster = cluster_eucs;
    kMeans.euc.aic = aic_euclid;
    kMeans.euc.bic = bic_euclid;
    kMeans.mah.group = group_mah;
    kMeans.mah.cluster = cluster_mahs;
    kMeans.mah.aic = aic_mah;
    kMeans.mah.bic = bic_mah;
    kMeans.mah.cov_f = cov_f;
    kMeans.type = {'Euclidean','Mahalanobis'};
end
kMeans.kmin = kmin;
kMeans.kmax = kmax;


% if (ma & eu)
%     save(resultfile, 'group_mah', 'group_euclid', 'cluster_mahs','cluster_eucs','aic_mah', 'bic_mah', ...
%         'aic_euclid', 'bic_euclid', 'cov_f', 'kmin', 'kmax');
% elseif (ma)
%     save(resultfile, 'group_mah','cluster_mahs', 'aic_mah', 'bic_mah', 'cov_f', 'kmin', 'kmax');
% else
%     save(resultfile, 'group_euclid','cluster_eucs', 'aic_euclid', 'bic_euclid', 'kmin', 'kmax');
% end