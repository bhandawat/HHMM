function [H,window] = hopkins(data,m,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hopkins Function for calculating the Hopkins statistics
% Implementation of the Hopkins statistics
%      (cf. "A New Method for determining the Type of Distribution of Plant
%       Individuals.", B. Hopkins and J. G. Skellam, 1954)
%
% Inputs:
%    data: nxd matrix where n = number of observations, d = number of
%    dimensions
%    m: number of points to use (m<<n)
%    type: 'Norm' = normal data
%          'Prob' = probability distribution

% Outputs:
%    H: hopkins statistic
%    window: hyperrectangle bounds (d x 2 matrix) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2018, Liangyu Tao

[n,dim] = size(data);

% use m<<n, if not user specified, use n = 2 order magnitude smaller
if isempty(m)
    m = ceil(n./100);
end

if strcmpi(type,'Prob')
    Y = drchrnd(ones(1,dim), m);
    window = ones(dim,2);
else
    % generate random points in a nDimensional window
    window = zeros(dim,2);
    Y = zeros(m,dim);
    for d = 1:dim
        window(d,:) = [min(data(:,d)) max(data(:,d))];
        center =  (max(data(:,d))+min(data(:,d)))./2;
        range =  (max(data(:,d))-min(data(:,d)))./2;

        Y(:,d) =  center+range.*(rand(m, 1)-0.5);
    end
end
X = data;

if strcmpi(type,'Prob')
    % use KL divergence if probability distribution
    % calculate the dist of each point in X from it's nearest neighbour in X
    d = zeros(length(X));
    for i = 1:length(X)
        for j = 1:length(X)
            d(i,j)=kullback_leibler_divergence(X(i,:),X(j,:));
        end
    end
    D = d+d';
    D(logical(eye(size(D)))) = 1;
    [wi,~] = min(D);
    
    % calculate the dist of each point in Y from it's nearest neighbour in X
    d = zeros(size(X,1),size(Y,1));
    d2 = zeros(size(X,1),size(Y,1));
    for i = 1:size(X,1)
        for j = 1:size(Y,1)
            d(i,j)=kullback_leibler_divergence(X(i,:),Y(j,:));
            d2(i,j)=kullback_leibler_divergence(Y(j,:),X(i,:));
        end
    end
    D2 = d+d2;
    D2(logical(eye(size(D2)))) = 1;
    [ui,~] = min(D2);
else
    % use euclidean distance for everything else
    % calculate the dist of each point in X from it's nearest neighbour in X
    D = squareform(pdist(X));
    D(logical(eye(size(D)))) = 1;
    [wi,~] = min(D);
    
    % calculate the dist of each point in Y from it's nearest neighbour in X
    D2 = pdist2(data,Y);
    [ui,~] = min(D2);
end


% sample m nearest neighbor pairs from wi;
sampNdx = sort(datasample([1:1:n],m));
wi = wi(sampNdx);

% calculate the Hopkins statistics
H = sum(ui)./(sum(ui)+sum(wi));

end