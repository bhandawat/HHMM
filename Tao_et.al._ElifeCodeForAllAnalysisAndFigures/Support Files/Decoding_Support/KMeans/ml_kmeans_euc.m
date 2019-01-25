function [cluster, group] = ml_kmeans_euc(feat, groupinfo , k, idx)
% FUNCTION [CLUSTER, GROUP] = ML_KMEANS_EUC(FEAT, GROUPINFO, K)
% K-means clustering using Euclidean distance and predefined group information
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
% Xiang Chen
% Aug 13, 2002

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

n = size(feat, 1);

if ((~isempty(groupinfo)) & (n ~= length(groupinfo)))
    error('Not equal size.');
end

%idx = randperm(n);
%idx = 1:n;
idx2(1) = idx(1);
for m = 2 : k
    dist = [];
    for o = 1 : length(idx2)
        dist(:, o) = sqrt(sum((feat - repmat(feat(idx2(o), :), n, 1)) .* (feat - repmat(feat(idx2(o), :), n, 1)), 2));
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
done = 0;
cluster = cell(1, k);
while (~done)
    %keyboard
    cluster = cell(1, k);
    %no_clust = zeros(1, k);
    for m = 1 : n
        index = 1;
        dmin = xc_squaresum2(feat(m, :) - center(1, :));
        for l = 2:k
            d = xc_squaresum2(feat(m, :) - center(l, :));
            if (d < dmin)
                %keyboard
                index = l;
                dmin = d;
            end
        end
        cluster{index} = [cluster{index} m];
        %keyboard
    end
    for l = 1 :k
        if (length(cluster{l} > 0))
            center(l, :) = mean(feat(cluster{l}, :), 1);
        end
    end
    sume = 0;
    for m = 1 : k
        for l = 1 : length(cluster{m})
            sume = sume + xc_squaresum2(feat(cluster{m}(l), :)...
                -center(m,:));
        end
    end
    if (~exist('prev', 'var'))
        prev = sume;
    else
        if (sume < prev)
            prev = sume;
        else
            done = 1;
        end
    end
end

if (~isempty(groupinfo))
    for m = 1 : length(cluster)
        group{m} = sort(groupinfo(cluster{m}));
    end
end

function x = xc_squaresum2(mtx)
% FUNCTION X = XC_SQUARESUM2(MTX)
% Calculate the suquaresum of each column of the vector.  If
% the input is a 2D matrix, each row will be calculated.

%%%%%%%commented by T. Zhao%%%%%%%%%
% x = zeros(size(mtx, 1), 1);
%
% for m = 1 : size(mtx, 1)
%      for n = 1 : size(mtx, 2)
%          x(m) = x(m) + mtx(m, n) * mtx(m, n);
%      end
%      %x(m) = sqrt(x(m));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%added by T. Zhao%%%%%%%%%%
x = sum(mtx.^2,2);