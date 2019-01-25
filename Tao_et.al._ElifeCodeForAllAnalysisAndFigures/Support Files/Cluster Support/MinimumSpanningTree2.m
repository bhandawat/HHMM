function [T_prime,rej0,setOne,setTwo,AllPts,mst] = MinimumSpanningTree2(data,dim,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinimumSpanningTree2 Function for conducting the MST based cluster
% analysis
% Implementation of the MST based cluster tendency
%      (cf. "Uniformity Testing Using Minimal Spanning Tree.",
%       A. Jain et al. , 2002)
%
% Inputs:
%    data: nxd matrix where n = number of observations, d = number of
%    dimensions (This is denoted as {x})
%    dim: number of dimensions to consider (dim<=d)
%    type: 'Norm' = normal data
%          'Prob' = probability distribution

% Outputs:
%    T_prime: z-score
%    rej0: percent of points that passes the edge inconsistancy test
%    setOne: ndx in AllPts for {x}
%    setTwo: ndx in AllPts for {y}
%    AllPts: (n+m)xdim matrix of all the data points used in this analysis
%    mts: minimum spanning tree

% *Note: Sets {x} and {y} are defined as half of set {data}. Therefore, {x}
% represents ~ half of the input data defined by a high average edge length
% and {y} represents ~ half of the input data defined by a low average edge
% length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2018, Liangyu Tao

if isempty(dim)
    dim = size(data,2);
end
dat = data(:,1:dim);

% number of input and test points
nData = size(data,1);
DataNdx = 1:nData;

% calculate the adjascency matrix and mst
[A,G,mst] = calcMST(dat,nData,dim,type);

% inconsistancy test, ro>6 to significantly reject the null (Ho)
lj = mst.Edges.Weight;
mu = mean(lj);
sigma = std(lj);
l_norm = (lj-mu)./sigma;
ro = max(l_norm);
rej0 = sum(l_norm>=6)./length(l_norm);

% degree of nodes
d = histcounts(mst.Edges.EndNodes,[1:1:nData+1]);

% compute high average edge length distribution (p_bar) and low average
% edge length distribution (bar_p)
wi_bar = zeros(1,nData);
bar_wi = zeros(1,nData);
for i = 1:nData
    [r,~] = find(mst.Edges.EndNodes==i);
    wi_bar(i) = sum(exp(l_norm(r)))./d(i);
    bar_wi(i) = sum(exp(-l_norm(r)))./d(i);
end
s_bar = sum(wi_bar);
bar_s = sum(bar_wi);

p_bar = wi_bar./s_bar;
bar_p = bar_wi./bar_s;

nr = floor(nData./2);
setOne = 1:1:nr;
setTwo = nr+1:1:2*nr;

% sample data without replacement based on p_bar and bar_p
x_bar = datasample(DataNdx,nr,'Weights',p_bar,'Replace',false);
bar_x = datasample(DataNdx,nr,'Weights',bar_p,'Replace',false);

% Run mst on the two subsamples
AllPts = [dat(x_bar,:);dat(bar_x,:)];

% calculate the adjascency matrix and mst
[A2,G2,mst2] = calcMST(AllPts,nr.*2,dim,type);

% check to see if the edges of the mts connect the same sets of
% points (xx,yy joints) or different sets (xy joints)
SameSetVal = sum(mst2.Edges.EndNodes<=nr,2);
xx_yy_joins = sum(SameSetVal~=1);
xy_joins = sum(SameSetVal==1);

% redefine variables
T = xy_joins;
L = nr*2;
m = nr;
n = nr;
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

function [A,G,mst] = calcMST(AllPts,nPts,dim,type)

% create a dissimilarity matrix from the entire set of points
if strcmpi(type,'Prob')
    % use KL divergence if probability distribution
    AllPts = AllPts+eps;
    AllPts = AllPts./repmat(sum(AllPts,2),1,dim);
    d = zeros(nPts);
    for i = 1:nPts
        for j = 1:nPts
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

end