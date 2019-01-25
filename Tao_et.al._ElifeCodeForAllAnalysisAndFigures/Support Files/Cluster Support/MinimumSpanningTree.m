function [T_prime,RandNdx,DataNdx,AllPts,K,inSamp,mst] = MinimumSpanningTree(data,dim,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinimumSpanningTree Function for conducting the MST based cluster
% analysis
% Implementation of the MST based cluster tendency
%      (cf. "Testing for uniformity in mutidimensional data.",
%       S. Smith and A. Jain, 1984)
%
% Inputs:
%    data: nxd matrix where n = number of observations, d = number of
%    dimensions (This is denoted as {x})
%    dim: number of dimensions to consider (dim<=d)
%    type: 'Norm' = normal data
%          'Prob' = probability distribution

% Outputs:
%    T_prime: z-score
%    RandNdx: ndx in AllPts for {y}
%    DataNdx: ndx in AllPts for {x}
%    AllPts: (n+m)xdim matrix of all the data points used in this analysis
%    K: jxdim matrix of indices of the points in {X} that generate the
%    convex hull
%    inSamp: mxd matrix of generated points in the convex hull
%    mts: minimum spanning tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2018, Liangyu Tao

if isempty(dim)
    dim = size(data,2);
end
dat = data(:,1:dim);

%--------------------------------------------------------------------------
if strcmpi(type,'Prob')
    % For probability distributions, sample from a Dirichlet distribution
    in = ones(1,dim);
    inSamp = drchrnd(ones(1,dim), size(data,1));
    K = [];
else
    % For normal data, find the convex hull over the data
    [K,~] = convhulln(dat,{'Qt','Qx'});
    
    % generate uniform, random points in a rectangular region and keep
    % the same amount of points (that lie in the convex hull) as the
    % original data
    X_cube = rand(3000, 3)-0.5;
    in = inhull(X_cube,dat,K,0);
    in = find(in,size(data,1));
    inSamp = X_cube(in,:);
end
%--------------------------------------------------------------------------

% number of input and test points
nRandPts = length(in);
nData = size(data,1);

RandNdx = 1:nRandPts;
DataNdx = nRandPts+1:nRandPts+nData;
AllPts = [inSamp;dat];

% create a dissimilarity matrix from the entire set of points
if strcmpi(type,'Prob')
    % use KL divergence if probability distribution
    AllPts = AllPts+eps;
    AllPts = AllPts./repmat(sum(AllPts,2),1,dim);
    d = zeros(length(AllPts));
    for i = 1:length(AllPts)
        for j = 1:length(AllPts)
            d(i,j)=kullback_leibler_divergence(AllPts(i,:),AllPts(j,:));
        end
    end
    A = d+d';
else
    % use euclidean distance for everything else
    D = pdist(AllPts);
    A = squareform(D);
end

% create a complete weighted graph of all the points
G = graph(A);
% calculate the minimum spanning tree
mst = minspantree(G);

% check to see if the edges of the mts connect the same sets of
% points (xx,yy joints) or different sets (xy joints)
SameSetVal = sum(mst.Edges.EndNodes<=nRandPts,2);
xx_yy_joins = sum(SameSetVal~=1);
xy_joins = sum(SameSetVal==1);

% redefine variables
T = xy_joins;
L = nRandPts+nData;
m = nRandPts;
n = nData;
% degree of nodes
d = histcounts(mst.Edges.EndNodes,[1:1:L+1]);

% calculate the expected mean and var of the number of xy joints
E_T = 2.*m.*n./L;
C = 0.5.*sum(d.*(d-1));
v_T = 2*m*n/(L*(L-1)) * ((2*m*n-L)/L +((C-L+2)/((L-2)*(L-3)))*(L*(L-1)-4*m*n+2));

% calculate the normalized version of the xy joints (~ follows standard
% normal distribution)
T_prime = (T-E_T)/sqrt(v_T);

end