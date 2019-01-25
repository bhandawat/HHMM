function [cluster, group] = ml_kmeans_mah(feat, groupinfo , k, f_cov)
% FUNCTION [CLUSTER, GROUP] = ML_KMEANS_MAH(FEAT, GROUPINFO, K, F_COV)
% K-means clustering using Mahalanobis distance and predefined group
% information
% (i.e. a few observations from the same cell clone and should belong to the
% same cluster)
% feat - a m*n feature matrix with m observations and n features
% groupinfo - a m-element array with group information, optional
% k - number of clusters
% cluster - a k-element cell array, each containing the indices of the 
%           observations belong to this cluster
% group - a k-element cell array similar to cluster, just replace the indices
%         to the corresponding information in groupinfo, only available when
%         groupinfo is supplied
% f_cov - an optional input of the covariance matrix.  It will be
%         calculated from the whole dataset if ignored
% Xiang Chen
% Xiang Chen Sept 25, 2003 (Modified from mv_mahalkmeans)
 
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

%No. of instances
n = size(feat, 1);

if ((~isempty(groupinfo)) & (n ~= length(groupinfo)))
     error('Not equal size.');
end

% Ranomize the initial seeds
idx = randperm(n);
idx2(1) = idx(1);
for m = 2 : k
    dist = [];
    if (exist('f_cov', 'var'))
      for n1 = 1 : length(idx2)
         for o = 1 : size(feat, 1)
             dist(o, n1) = (feat(o, :) - feat(idx2(n1))) * inv(f_cov) *(feat(o, :) - feat(idx2(n1)))';
         end
      end
    else 
        dist = squareform(pdist(feat, 'mahal'));
        dist = dist(:, idx2);
    end
    [v, I] = sort(sum(dist, 2));
    done = 0;
    c = size(dist, 1);
    while (~done)
        if (find(idx2 == I(c)))
            c = c - 1;
        else
            idx2(m) = I(c);
            done = 1;
        end
    end
end
center = feat(idx2, :);

clust = zeros(n, k);
id = eye(k);
clust(idx2,:) = id;
done = 0;
while (~done)
	dist = [feat
            center];
    cen2 = center;
    if (exist('f_cov', 'var'))
       for m = 1 : k
         for o = 1 : size(feat, 1)
             pd(m, o) = (feat(o, :) - center(m, :)) * inv(f_cov) * (feat(o, :) - center(m, :))';
         end
      end
    else 
        pd = squareform(pdist(dist, 'mahal'));
        pd = pd(n + 1 : n + k, 1 : n);
    end
    [Y, I] = min(pd, [], 1);
    clust = id(I, :);
    no_instances = sum(clust, 1);

    for m = 1 : k
	if (no_instances(m) > 0)
	    center(m, :) = mean(feat(find(clust(:,m)),:), 1);
        end
    end
	
    if (center == cen2)
	done = 1;
    end
end

for m = 1 : k
	tmp = clust(:, m);
	cluster{m} = find(tmp);
end

 
if (~isempty(groupinfo))
     for m = 1 : length(cluster)
         group{m} = sort(groupinfo(cluster{m}));
     end
end