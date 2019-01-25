function [data,model,likely_high_state_sorted,kMNew] = GLMClusterAnalysis(kMeans,data,model,likely_high_state_sorted,c2cons)
%GLMClusterAnalysis Takes in kMeans clusters and uses GLM to redistribute
%the three flies that most likely don't belong to the two most represented
%clusters. This finds an improved local minimum over kMeans

% Inputs:
%    kMeans: Cell containing kMeans clustering of flies.
%    data: Data cell array containing raw observables
%    model: HHMM information
%    likely_high_state_sorted: nx10799 matrix of sorted HLS tracks
%
% 2017, Liangyu Tao

delete(gcp('nocreate'))
parpool('local',2)

fit = cell(2,2);
HLSPredict = cell(2,2);
Ctemp = cell(2,5);
CMissingTemp = cell(2,5);
cNum = find(~cellfun(@isempty,c2cons));

% Inside
loc = 1;
[Ctemp,CMissingTemp,fit,HLSPredict,clusts{loc}] = calcGLMFit(data,model,...
    likely_high_state_sorted,kMeans,Ctemp,CMissingTemp,fit,HLSPredict,loc,c2cons);

% Outside
loc = 2;
[Ctemp,CMissingTemp,fit,HLSPredict,clusts{loc}] = calcGLMFit(data,model,...
    likely_high_state_sorted,kMeans,Ctemp,CMissingTemp,fit,HLSPredict,loc,c2cons);
%load('tmp.mat')

%fit2 = fit;
kM = kMeans;
ndxBF = zeros(2,2);
HLSPredictFit = [];
for loc = 1:2
    n = find(~cellfun(@isempty,Ctemp(loc,:)));
    for c = 1:length(n)
        [~,ndx] = sort(fit{loc,c});
        fit{loc,c} = fit{loc,c}(ndx);
        HLSPredict{loc,c} = HLSPredict{loc,c}(ndx,:);
        Ctemp{loc,n(c)} = Ctemp{loc,n(c)}(ndx,:);
        CMissingTemp{loc,n(c)} = CMissingTemp{loc,n(c)}(ndx,:);
        
        ndxBF(loc,c) = find(fit{loc,c}(:,1)==max(fit{loc,c}(:,1)));
        HLSPredictTmp = HLSPredict{loc,c}(ndxBF(loc,c),:);
        HLSPredictFit{loc,c} = HLSPredictTmp;
        
        for i = 1:1:length(ndx)
            HLSPredictTmp = HLSPredict{loc,c}(ndx(i),:);
            if min(HLSPredictTmp)>0.5
                ndxBF(loc,c) = i;
                HLSPredictFit{loc,c} = HLSPredictTmp;
            end
        end
    end
end

loc = 1;clust_to_cons = clusts{loc};
numCluster = kMeans.data{loc,3}.optCluster;
CMissing = cell(1,2);
for c = 1:length(clust_to_cons)
    kM.data{loc,3}.cluster{numCluster}{clust_to_cons(c)}=Ctemp{loc,clust_to_cons(c)}(ndxBF(loc,c),:);
    CMissing{loc} = [CMissing{loc},CMissingTemp{loc,clust_to_cons(c)}(ndxBF(loc,c),:)];
end

loc = 2;clust_to_cons = clusts{loc};
numCluster = kMeans.data{loc,3}.optCluster;
for c = 1:length(clust_to_cons)
    kM.data{loc,3}.cluster{numCluster}{clust_to_cons(c)}=Ctemp{loc,clust_to_cons(c)}(ndxBF(loc,c),:);
    CMissing{loc} = [CMissing{loc},CMissingTemp{loc,clust_to_cons(c)}(ndxBF(loc,c),:)];
end


for loc =1:2
    numCluster = kM.data{loc,3}.optCluster;
    currentCluster = kM.data{loc,3}.cluster{numCluster};
    
    possClust = repmat({1:numCluster},1,6);
    comb = cell(1, numel(possClust)); %set up the varargout result
    [comb{:}] = ndgrid(possClust{:});
    comb = cellfun(@(x) x(:), comb,'uniformoutput',false); %there may be a better way to do this
    Perm = [comb{:}];
    
    N = size(Perm,1);
    ppm = ParforProgMon('Processing: ', N , 1);
    parfor c = 1:size(Perm,1)
        kMTemp = kM;
        for c2 = unique(Perm(c,:))
            kMTemp.data{loc,3}.cluster{numCluster}{c2} = [kMTemp.data{loc,3}.cluster{numCluster}{c2} CMissing{loc}(Perm(c,:)==c2)];
        end
        [~,~,~,PCHLSALL,~,~] = GLMCluster(data,model,likely_high_state_sorted,kMTemp,c2cons);
        tmp2(c,:) = PCHLSALL{cNum,loc}.KC(loc,:);
        ppm.increment();
    end
    HLSPredict{loc,1} = tmp2;
    AllPerm{loc} = Perm;
    clearvars tmp2
end
%load('temp2.mat')

kMNew = kM;
HLSPredict2 = [];
for loc = 1:2
    numCluster = kM.data{loc,3}.optCluster;
    for i = 1:size(AllPerm{loc},1)
        HLSPredict2{loc,1}(i,:) = cell2mat(HLSPredict{loc,1}(i,:));
    end
    goodComb = find(min(HLSPredict2{loc}')'>0.5);
    if isempty(goodComb)
        goodComb = find(min(HLSPredict2{loc}')'==max(min(HLSPredict2{loc}')'));
    end
    fit2Global = sum(HLSPredict2{loc}(goodComb,:),2);
    likelyClusterNdx = goodComb(fit2Global==max(fit2Global));
    likelyCluster = AllPerm{loc}(likelyClusterNdx,:);
    
    for c = unique(likelyCluster)
        kMNew.data{loc,3}.cluster{numCluster}{c} = [kMNew.data{loc,3}.cluster{numCluster}{c} CMissing{loc}(likelyCluster==c)];
    end
end
end

function [Ctemp,CMissingTemp,fit,HLSPredict,clust_to_cons] = calcGLMFit(data,model,...
    likely_high_state_sorted,kMeans,Ctemp,CMissingTemp,fit,HLSPredict,loc,c2cons)

numCluster = kMeans.data{loc,3}.optCluster;
currentCluster = kMeans.data{loc,3}.cluster{numCluster};
clust_to_cons = find(cellfun(@length,currentCluster)>=10);
nPerClust = floor(6./length(clust_to_cons));
cNum = find(~cellfun(@isempty,c2cons));
tic
for i = 1:length(clust_to_cons)
    [~,~,~,PCHLSALL,~,~] = GLMCluster(data,model,likely_high_state_sorted,kMeans,c2cons);
    cClust = clust_to_cons(i);
    [~,I]=sort(PCHLSALL{cNum,loc}.KC{loc,cClust},'ascend');
    I = currentCluster{cClust}(I);
    C = combnk(I(1:min(10,length(currentCluster{cClust}))),min(length(currentCluster{cClust}),10)-nPerClust);
    C = [C repmat(I(11:end),size(C,1),1)];
    
    for j = 1:size(C,1)
        CMissingTemp{loc,cClust}(j,:) = setdiff(I,C(j,:));
    end
    Ctemp{loc,cClust} = C;
    tmp = zeros(size(C,1),1);tmp2 = zeros(size(C));
    
    N = size(C,1);
    ppm = ParforProgMon('Processing: ', N , 1);
    parfor c = 1:size(C,1)
        kM = kMeans;
        kM.data{loc,3}.cluster{numCluster} = num2cell(ones(1,numCluster));
        kM.data{loc,3}.cluster{numCluster}{cClust} = C(c,:);
        [~,~,~,PCHLSALL,~,~] = GLMCluster(data,model,likely_high_state_sorted,kM,c2cons);
        tmp(c,1) = sum(PCHLSALL{cNum,loc}.KC{loc,cClust});
        tmp2(c,:) = PCHLSALL{cNum,loc}.KC{loc,cClust};
        ppm.increment();
    end
    fit{loc,i} = tmp;
    HLSPredict{loc,i} = tmp2;
end
toc

end
