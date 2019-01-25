function [] = subSamplingHLSDist(data_time,model,likely_high_state_sorted,c2cons,options,tBins,cNum,fig_title)

KMDatFile = options.empKMFile;
SynDatFile = options.synthFile;

if isempty(cNum)
    cNum = 1;
end
if isempty(tBins)
    tBins = [5 30:30:5400];
end

load(KMDatFile,'kMeans');
load(SynDatFile,'B_new');
B_new = B_new{cNum};
kMeans = kMeans{cNum};

nFlys = length(c2cons{cNum});
nIt = 30;
nHL = model.Qdim;
sce2Cons = {'b_o','d_i'};
[clustCent,clabels] = getClustCent(B_new,kMeans,nFlys);

% perform subsampling
d = [];
for t = 1:length(tBins)
    distRAll = [];distRBlock = [];distFirst = [];distLast = [];
    for fly = 1:nFlys
        idxoff=find(data_time{c2cons{cNum}(fly)}(11,1:end-1)<1 & ...
            data_time{c2cons{cNum}(fly)}(12,1:end-1)<1);                    % odor off and fly outside
        idxon=find(data_time{c2cons{cNum}(fly)}(11,1:end-1)==1 & ...
            data_time{c2cons{cNum}(fly)}(12,1:end-1)==1);                   % odor on and fly inside
        tmpHLS = likely_high_state_sorted(c2cons{cNum}(fly),:);
        HLSOff = tmpHLS(idxoff);
        HLSOn = tmpHLS(idxon);
        
        X = diff(HLSOff)~=0;
        BOff = find([true,X]);
        EOff = find([X,true]);
        DOff = EOff-BOff+1;
        
        X = diff(HLSOn)~=0;
        BOn = find([true,X]);
        EOn = find([X,true]);
        DOn = EOn-BOn+1;
        
        % for randomly generated indexes
        for n = 1:nIt
            if tBins(t)>=length(idxoff)
                distRAll{1}(:,n) = histcounts(tmpHLS(idxoff),[1:1:nHL+1],'Normalization','Probability');
            else
                tmp = [];
                while length(tmp)<tBins(t)
                    tmpNdx = datasample(BOff,1);
                    tmp = [tmp HLSOff(tmpNdx).*ones(1,DOff(BOff==tmpNdx))];
                end
                distRAll{1}(:,n) = histcounts(tmp(1:tBins(t)),[1:1:nHL+1],'Normalization','Probability');
            end
            if tBins(t)>=length(idxon)
                distRAll{2}(:,n) = histcounts(tmpHLS(idxon),[1:1:nHL+1],'Normalization','Probability');
            else
                tmp = [];
                while length(tmp)<tBins(t)
                    tmpNdx = datasample(BOn,1);
                    tmp = [tmp HLSOn(tmpNdx).*ones(1,DOn(BOn==tmpNdx))];
                end
                distRAll{2}(:,n) = histcounts(tmp(1:tBins(t)),[1:1:nHL+1],'Normalization','Probability');
            end
        end
        
        % for random blocks of bins
        if length(idxoff)-tBins(t)>1
            y = datasample(1:length(idxoff)-tBins(t),nIt);
            y2 = mat2cell(tmpHLS(y'+[1:1:tBins(t)]),ones(1,nIt),tBins(t));
        else
            y2 = mat2cell(repmat(tmpHLS(idxoff),nIt,1),ones(1,nIt),length(idxoff));
        end
        [hcell,~] = cellfun(@(x) histcounts(x,[1:1:nHL+1],'Normalization','Probability'), y2, 'Uni',0);
        distRBlock{1} = cell2mat(hcell)';
        
        if length(idxon)-tBins(t)>1
            y = datasample(1:length(idxon)-tBins(t),nIt);
            y2 = mat2cell(tmpHLS(y'+[1:1:tBins(t)]),ones(1,nIt),tBins(t));
        else
            y2 = mat2cell(repmat(tmpHLS(idxon),nIt,1),ones(1,nIt),length(idxon));
        end
        [hcell,~] = cellfun(@(x) histcounts(x,[1:1:nHL+1],'Normalization','Probability'), y2, 'Uni',0);
        distRBlock{2} = cell2mat(hcell)';
        
        
        % for the first n frames and the last n frames
        y = idxoff(1:min(tBins(t),length(idxoff)));
        distFirst{1}(:,1) = histcounts(tmpHLS(y),[1:1:nHL+1],'Normalization','Probability');
        y = idxoff(max(length(idxoff)-tBins(t),1):end);
        distLast{1}(:,1) = histcounts(tmpHLS(y),[1:1:nHL+1],'Normalization','Probability');
        
        
        y = idxon(1:min(tBins(t),length(idxon)));
        distFirst{2}(:,1) = histcounts(tmpHLS(y),[1:1:nHL+1],'Normalization','Probability');
        y = idxon(max(length(idxon)-tBins(t),1):end);
        distLast{2}(:,1) = histcounts(tmpHLS(y),[1:1:nHL+1],'Normalization','Probability');
        
        
        for sce = 1:2
            dTmp = zeros(length(clustCent.(sce2Cons{sce})),nIt);dTmp2 = zeros(length(clustCent.(sce2Cons{sce})),nIt);
            dTmp3 = zeros(length(clustCent.(sce2Cons{sce})),1);dTmp4 = zeros(length(clustCent.(sce2Cons{sce})),1);
            for c = 1:length(clustCent.(sce2Cons{sce}))
                for n = 1:nIt
                    dTmp(c,n) = sqrt(sum((clustCent.(sce2Cons{sce}){c}-distRAll{sce}(:,n)).^2));
                    dTmp2(c,n) = sqrt(sum((clustCent.(sce2Cons{sce}){c}-distRBlock{sce}(:,n)).^2));
                end
                
                dTmp3(c) = sqrt(sum((clustCent.(sce2Cons{sce}){c}-distFirst{sce}).^2));
                dTmp4(c) = sqrt(sum((clustCent.(sce2Cons{sce}){c}-distLast{sce}).^2));
                
            end
            [~,cAssign] = min(dTmp);
            [~,cAssign2] = min(dTmp2);
            [~,cAssign3] = min(dTmp3);
            [~,cAssign4] = min(dTmp4);
            
            d.RAll{sce,t}(fly,:) = cAssign==clabels.(sce2Cons{sce})(fly);
            d.RBlock{sce,t}(fly,:) = cAssign2==clabels.(sce2Cons{sce})(fly);
            d.First{sce,t}(fly,:) = cAssign3==clabels.(sce2Cons{sce})(fly);
            d.Last{sce,t}(fly,:) = cAssign4==clabels.(sce2Cons{sce})(fly);
        end
    end
    
end

RAll = zeros(2,length(tBins));RBlock = zeros(2,length(tBins));
First = zeros(2,length(tBins));Last = zeros(2,length(tBins));
for sce = 1:2
    for t = 1:length(tBins)
        RAll(sce,t) = mean(sum(d.RAll{sce,t}));
        RBlock(sce,t) = mean(sum(d.RBlock{sce,t}));
        First(sce,t) = mean(sum(d.First{sce,t}));
        Last(sce,t) = mean(sum(d.Last{sce,t}));
    end
end

% approximate chance from area of cluster and number of flies within it
features = drchrnd(ones(1,model.Qdim), 100000);
nSamp = size(features,1);nSce = length(sce2Cons);

pFly = zeros(nSce,5);pArea = zeros(nSce,5);
for i = 1:nSce
    cCent = clustCent.(sce2Cons{i});
    cCent2 = permute(cell2mat(cCent)',[3 2 1]);
    nC = size(cCent2,3);
    
    features2 = repmat(features,1,1,nC);
    
    tmp = bsxfun(@minus,features2,cCent2);
    dist = squeeze(cellfun(@norm,num2cell(tmp,2)));
    [~,N] = min(dist,[],2);
    
    for c = 1:nC
        pFly(i,c) = sum(clabels.(sce2Cons{i})==c)./nFlys;
        pArea(i,c) = sum(N==c)./nSamp;
    end
end

% approximate chance
pChance = sum(pFly.*pArea,2);

% plotting figure 7-S4
figure;set(gcf,'position',[855 49 824 918])
subplot(3,1,1);
plot(tBins./30,RAll'./nFlys);hold on;
plot([tBins(1)./30 tBins(end)./30],[pChance(1) pChance(1)],'--k')
plot([tBins(1)./30 tBins(end)./30],[pChance(2) pChance(2)],'--r')
ylim([0 1])
xlabel('time (s)');ylabel('probability correctly labeled');
title('Random Bins')
legend({'Before Outside','During Inside','Chance Outside','Chance Inside'},'Location','SouthEast')
subplot(3,1,2);
plot(tBins./30,RBlock'./nFlys);hold on
plot([tBins(1)./30 tBins(end)./30],[pChance(1) pChance(1)],'--k')
plot([tBins(1)./30 tBins(end)./30],[pChance(2) pChance(2)],'--r')
ylim([0 1])
xlabel('time (s)');ylabel('probability correctly labeled');
title('Random Blocks')
legend({'Before Outside','During Inside','Chance Outside','Chance Inside'},'Location','SouthEast')
subplot(3,1,3);
plot(tBins./30,[First;Last]'./nFlys);hold on
plot([tBins(1)./30 tBins(end)./30],[pChance(1) pChance(1)],'--k')
plot([tBins(1)./30 tBins(end)./30],[pChance(2) pChance(2)],'--r')
ylim([0 1])
xlabel('time (s)');ylabel('probability correctly labeled');
title('Random Blocks')
legend({'Before Outside Init','During Inside Init','Before Outside Last','During Inside Last','Chance Outside','Chance Inside'},'Location','SouthEast')

if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

end

function [sce,clabels] = getClustCent(B_new,kMeans,nFlys)
sce = [];clabels = [];
for state = 1:2
    for loc = 1:2
        clabels.(kMeans.scenario{loc,state}) = zeros(nFlys,1);
        optCluster=kMeans.maxCluster.data{loc,state}.optCluster;
        for i = 1:optCluster
            cluster=kMeans.maxCluster.data{loc,state}.cluster{1,optCluster}{1,i};
            clabels.(kMeans.scenario{loc,state})(cluster) = i;
            raw = B_new{loc,state}(:,cluster);
            sce.(kMeans.scenario{loc,state}){i} = mean(raw,2);
        end
    end
end
end