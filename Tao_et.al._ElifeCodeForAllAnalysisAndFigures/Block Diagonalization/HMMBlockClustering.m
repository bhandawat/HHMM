function [dur_HMM_sorted] = HMMBlockClustering(model,data,HMMsortNdx,likely_state_by_fly_sorted,param,fig_title)
K = param.K;
maxIt = param.maxIt;
tau = param.tau;
beta = param.beta;

TP = model.HMMs{1}.A.alpha-model.HMMs{1}.A.alpha_0;
TP_sorted=TP(HMMsortNdx,HMMsortNdx);
nGood = sum(sum(TP)>0);
TP_good = TP_sorted(1:nGood,1:nGood);
TP_good(TP_good<1) = 0;

TP_tot = repmat(sum(TP_good,1),nGood,1);
TP_normalized = (TP_good./TP_tot)';

figTitle = ['Clustering_tau' num2str(tau) '_beta' num2str(beta)];

[~,TP_sorted,~,~,HLSNdx] = PIB(TP_normalized,K,maxIt,tau,beta);

if length(HLSNdx)==K
figure;imagesc(TP_sorted);hold on
x = 0.5;y = 0.5;
for clust = 1:length(HLSNdx)
    sz = length(HLSNdx{clust});
    rectangle('Position',[x y sz sz],'LineWidth',2,'EdgeColor','w')
    x = x+sz;y = y+sz;
end
title('HMM Transition Probability w/ Clusters')
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

HLTP = zeros(10);
for clust = 1:length(HLSNdx)
    for clust2 = 1:length(HLSNdx)
        HLTP(clust,clust2) = sum(sum(TP_good(HLSNdx{clust},HLSNdx{clust2})));
    end
end

TPHL_tot = repmat(sum(HLTP,1),length(HLSNdx),1);
TPHL_normalized = HLTP./TPHL_tot;
TPHL_normalized = TPHL_normalized;

figure;imagesc(TPHL_normalized);title('Cluster Transition Probability')
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
likelyHLS = zeros(size(likely_state_by_fly_sorted));
for k = 1:length(HLSNdx)
    for n = 1:length(HLSNdx{k})
        likelyHLS(likely_state_by_fly_sorted==HLSNdx{k}(n)) = k;
    end
end

dur = cell(model.dim,1);
Track = cell(model.dim,34);
for fly = 1:34
    d1 = data{fly};
    for state = 1:10
        temp = likelyHLS(fly,:)==state;
        
        startNdx = find(diff([false temp])==1);
        endNdx = find(diff([temp false])==-1);
        
        if length(startNdx)>length(endNdx)
            endNdx = [endNdx 10799];
        end
        
        if length(endNdx)>length(startNdx)
            startNdx = [1 startNdx];
        end
        
        dur{state} = [dur{state} endNdx-startNdx+1];
        meanspeed = [];meancurvature = [];
        for t = 1:length(startNdx)
            Track{state,fly}{t} = d1(:,startNdx(t):endNdx(t));
            if endNdx(t)-startNdx(t)>10
                currenttrack = d1(:,startNdx(t):endNdx(t));
                longTrack{state,fly}{t} = currenttrack;
                
                meanspeed=[meanspeed nanmean(currenttrack(7,:))];
                meancurvature=[meancurvature nanmean(currenttrack(8,:))];
            end
        end
        avgspeed{fly,state}=meanspeed';
        avgcurvature{fly,state}=meancurvature';
    end
end
% Define the ratio as speed over curvature
for state = 1:10
    ratio(state) = mean(cell2mat(avgspeed(:,state)))./std(cell2mat(avgcurvature(:,state)));
end
[~,ndx] = sort(ratio);
HLSNdx = HLSNdx(ndx);

longTrack = longTrack(ndx,:);
dur = dur(ndx,:);

figure;imagesc(TP_sorted);hold on
x = 0.5;y = 0.5;
for clust = 1:length(HLSNdx)
    sz = length(HLSNdx{clust});
    rectangle('Position',[x y sz sz],'LineWidth',2,'EdgeColor','w')
    x = x+sz;y = y+sz;
end
title('HMM Transition Probability w/ Clusters')
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

HLTP = zeros(10);
for clust = 1:length(HLSNdx)
    for clust2 = 1:length(HLSNdx)
        HLTP(clust,clust2) = sum(sum(TP_good(HLSNdx{clust},HLSNdx{clust2})));
    end
end

TPHL_tot = repmat(sum(HLTP,1),length(HLSNdx),1);
TPHL_normalized = HLTP./TPHL_tot;
TPHL_normalized = TPHL_normalized;

figure;imagesc(TPHL_normalized);title('Cluster Transition Probability')
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% calculate duration
durHist = zeros(10,2000);
maxTrack = zeros(10,1);
for state = 1:10
    temp = dur{state};
    if ~isempty(temp)
        n = histcounts(temp,max(temp));
        durHist(state,1:max(temp)) = n;
        maxTrack(state) = max(temp);
    end
end

dur_HMM_sorted = dur(1:10);

forwardTracks(longTrack,maxTrack)
suptitle('Cluster Forward Tracks')
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end
HLSObsModContour2(data,model,likelyHLS,HMMsortNdx,HLSNdx)
suptitle('Cluster Contours')
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

save('mat_Files/HHMM_HMM_comparison.mat','data','model','TP_sorted','HLSNdx','longTrack','maxTrack','likelyHLS','HMMsortNdx','HLSNdx','param')
end
end
