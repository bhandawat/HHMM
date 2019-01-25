function [NumPoints,durVec,indDurVec,type,HLSDistRawOn,HLSDistRawOff,PCHLSComp,...
    bs_all,bad_Fly_List,flyClustNdx] = GLM_decode(data,model,likely_high_state_sorted,c2cons)
%GLM_decode Takes in the HLS for each fly track to fit to a logistic
%regression model as a means of decoding Fly behavior.

% Inputs:
%    data: 1 x nFly cell array with empirical fly data
%    model: HHMM model
%    likely_high_state_by_fly: 34x10799 matrix of unsorted HLS tracks

% Outputs:
%    NumPoints: 1x2 cell array (1 = on, 2 = off) for avg number of points
%    in scenario. (in each cell, c1 = in, c2 = out)
%    durVec: 2x2 cell array of duration of contiguous segments of HLS
%    indDurVec: Same as durVec, but in each cell array is a nHLS x nFly
%    cell array
%    type: matrix containing GLM results
%    HLSDistRawOn: empirical HLS distributions for during odor
%    HLSDistRawOff: empirical HLS distributions for before odor
%    PCHLSComp: HLS distribution matrix for HLS dist after PCA
%    bs_all: array showing number of grouped points used for GLM
%    bad_Fly_List: array showing number of points considered for GLM
%    flyClustNdx: 1xnCluster cell array with fly index

% Note that the output of this will go into other analysis
%
% 2017, Bhandawat Lab

close all
dc=0:0.01:1;
[~,idx]=max(model.p);
plotidx=0;
opts=statset('glmfit');
opts.MaxIter=1000;
% number of points to average over
bs_all = 30;
data_old = data;
data = post_hoc_preprocess_LT(data_old);
warning('off')
lastwarn('')
warnArray = [];
warnArray.group = cell(1,35);
warnArray.ind = cell(1,34);
bad_points = [];
bad_points.j = [];
bad_points.l = [];
NumPoints = cell(1,2);
flyClustNdx = cell(1,model.NC);


durVec = cell(2,2);
indDurVec = cell(2,2);
for l = 1:4
    durVec{l}=cell(1,model.Qdim);
    indDurVec{l}=cell(model.Qdim,34);
end

for m = 1:length(bs_all)
    bs = bs_all(m);
    for c=1:length(c2cons)
        i = c2cons(c);
        for l = 1:2
            X{l}{i}=[];                                                     %
            Y{l}{i}=[];                                                     %
            rawX{l}{i}=[];                                                  % Speed and Curvature Predictors
            rawY{l}{i}=[];                                                  % Speed and Curvature Descriptor
            BX{l}{i}=[];                                                    % High Level States Predictors
            BY{l}{i}=[];                                                    % High Level States Descriptor
            TPX{l}{i}=[];                                                   %
            TPY{l}{i}=[];                                                   %
            idxoff{l} = [];                                                 % Index for stim off                                                %
            idxon{l} = [];                                                  % Index for stim on
            BXHLSComp{i}{m}{l} = [];
            TPXTPComp{i}{m}{l} = [];
        end
        idx2=find(idx==i);  % idx of datapoints assigned to cluster i
        flyClustNdx{i} = idx2;
        np=length(idx2);
        if(np>0)
            plotidx=plotidx+1;
            bad_Fly_List = [];
            px=ceil(sqrt(np));
            py=ceil(np/px);
            for j=1:np
                flyN = idx2(j);
                d1=data{flyN};
                o1=d1(12,:);
                
                p1=model.HHMMs{i}.p{flyN};                               % LL state probability
                clear p2
                for k=1:model.Qdim
                    p2unsorted(k,:)=sum(p1(model.HHMMs{i}.Aidx{k},:),1);            % HL state probability
                end
                
                p2 = zeros(10,10799);
                for k=1:model.Qdim
                    ndx = find(likely_high_state_sorted(flyN,:)==k);
                    p2(k,ndx) = 1;
                end
                
                TP = model.HHMMs{i}.getQxi(data(flyN));
                TP = TP{1};
                TP = reshape(TP,model.Qdim^2,size(TP,3));
                
                idxoff{2}=find(d1(11,1:end-1)<1 & d1(12,1:end-1)<1);        % odor off and fly outside
                idxoff{1}=find(d1(11,1:end-1)<1 & d1(12,1:end-1)==1);       % odor off and fly inside
                idxon{2}=find(d1(11,1:end-1)==1 & d1(12,1:end-1)<1);        % odor on and fly outside
                idxon{1}=find(d1(11,1:end-1)==1 & d1(12,1:end-1)==1);       % odor on and fly inside
                
                for l = 1:2
                    TransitionIdxOn = find(diff(idxon{l})>1);
                    TransitionIdxOff = find(diff(idxoff{l})>1);
                    
                    ponAll=p2(:,idxon{l});                                  % Prob of HL states given stim on
                    poffAll=p2(:,idxoff{l});
                    ponVec = likely_high_state_sorted(flyN,idxon{l});         % ponAll with same formate as likely_high_state_sorted
                    poffVec = likely_high_state_sorted(flyN,idxoff{l});
                    TPonAll = TP(:,idxon{l});
                    TPoffAll = TP(:,idxoff{l});
                    ndx2 = sum(ponAll,1)==0;
                    ndx22 = find(ndx2 == 1);
                    for tp = 1:length(TransitionIdxOn)
                        TransitionIdxOn(tp) = TransitionIdxOn(tp)-length(find(ndx22<=TransitionIdxOn(tp)));
                    end
                    
                    ponAll(:,ndx2) = [];
                    ponVec(:,ndx2) = [];
                    TPonAll(:,ndx2) = [];
                    
                    ndx2 = sum(poffAll,1)==0;
                    ndx22 = find(ndx2 == 1);
                    for tp = 1:length(TransitionIdxOff)
                        TransitionIdxOff(tp) = TransitionIdxOff(tp)-length(find(ndx22<=TransitionIdxOff(tp)));
                    end
                    
                    poffAll(:,ndx2) = [];
                    poffVec(:,ndx2) = [];
                    TPoffAll(:,ndx2) = [];
                    NumPoints{1}(flyN,l) = size(ponAll(:,1:end),2);
                    NumPoints{2}(flyN,l) = size(poffAll(:,1:end),2);
                    HLSDistRawOn{flyN,l} = mean(ponAll(:,1:end),2);
                    HLSDistRawOff{flyN,l} = mean(poffAll(:,1:end),2);
                    mean(ponAll(:,1:end),2);
                    
                    TransitionIdxOn = [0 TransitionIdxOn length(ponVec)];
                    TransitionIdxOn = unique(TransitionIdxOn);
                    for tp = 1:length(TransitionIdxOn)-1
                        track = ponVec(TransitionIdxOn(tp)+1:TransitionIdxOn(tp+1));
                        startIndex = [1, find(diff(track)~=0)+1];
                        startIndex = unique(startIndex);
                        endIndex = [startIndex(2:end) length(track)];
                        for iiii = 1:length(startIndex)
                            initState = track(startIndex(iiii));
                            durVec{1,l}{initState} = [durVec{1,l}{initState} endIndex(iiii)-startIndex(iiii)];
                            indDurVec{1,l}{initState,flyN} = [indDurVec{1,l}{initState,flyN} endIndex(iiii)-startIndex(iiii)];
                        end
                    end
                    
                    TransitionIdxOff = [0 TransitionIdxOff length(poffVec)];
                    TransitionIdxOff = unique(TransitionIdxOff);
                    for tp = 1:length(TransitionIdxOff)-1
                        track = poffVec(TransitionIdxOff(tp)+1:TransitionIdxOff(tp+1));
                        startIndex = [1, find(diff(track)~=0)+1];
                        startIndex = unique(startIndex);
                        endIndex = [startIndex(2:end)-1 length(track)];
                        endIndex = unique(endIndex);
                        for iiii = 1:length(startIndex)
                            initState = track(startIndex(iiii));
                            durVec{2,l}{initState} = [durVec{1,l}{initState} endIndex(iiii)-startIndex(iiii)];
                            indDurVec{2,l}{initState,flyN} = [indDurVec{2,l}{initState,flyN} endIndex(iiii)-startIndex(iiii)];
                        end
                    end
                    
                    for kk = 1:model.Qdim
                        tmp = TPonAll((kk-1)*model.Qdim+1:kk*model.Qdim,:);
                        tmp2 = TPoffAll((kk-1)*model.Qdim+1:kk*model.Qdim,:);
                        HLSDistTPon{flyN,l,kk} = sum(tmp,2)./sum(tmp(:));
                        HLSDistTPoff{flyN,l,kk} = sum(tmp2,2)./sum(tmp2(:));
                    end
                end
                
                %if (sum(cellfun(@isempty,idxon)) + sum(cellfun(@isempty,idxoff)))>0
                %remove flies with not enough points (361 > 4*90 ==> at least 4 points in average over 90)
                if length(idxon{1})<91 || length(idxoff{1})<91 || length(idxon{2})<91 || length(idxoff{2})<91
                    display(['Fly ' num2str(flyN) ' is a bad fly'])
                    display([num2str(length(idxon{1})) '  ' num2str(length(idxoff{1}))...
                        '  ' num2str(length(idxon{2})) '  ' num2str(length(idxoff{2}))])
                    bad_Fly_List = [bad_Fly_List, flyN];
                else
                    for l = 1:2                                             %l = 1 for inside, l = 2 for outside odor ring
                        
                        %setting chance to zero by tossing out points until
                        %idxon and idxoff are same length
                        pon=p2(:,idxon{l});                                 % Prob of HL states given stim on
                        poff=p2(:,idxoff{l});                               % Prob of HL states given stim off
                        ndx2 = sum(pon,1)==0;
                        idxon{l}(ndx2) = [];
                        pon(:,ndx2) = [];
                        ndx2 = sum(poff,1)==0;
                        idxoff{l}(ndx2) = [];
                        poff(:,ndx2) = [];
                        len=min(size(pon,2),size(poff,2));
                        pon=pon(:,1:len);
                        poff=poff(:,1:len);
                        
                        %len=min(length(idxon{l}),length(idxoff{l}));
                        idxoff{l}=idxoff{l}(1:len);
                        idxon{l}=idxon{l}(1:len);
                        
                        d2on=d1(7:8,idxon{l})';
                        d2on(8,:)=abs(d2on(8,:));                           % speed and curvature given stim on
                        d2off=d1(7:8,idxoff{l})';
                        d2off(8,:)=abs(d2off(8,:));                         % speed and curvature given stim off
                        
                        % this section is added because some data points were counted multiple times
                        pon2 = zeros(model.Qdim, floor(size(pon,2)/bs));
                        poff2 = zeros(model.Qdim, floor(size(pon,2)/bs));
                        d22on = zeros(floor(size(pon,2)/bs),2);
                        d22off = zeros(floor(size(pon,2)/bs),2);
                        % no overlap
                        for k=1:floor(size(pon,2)/bs)
                            pon2(:,k)=mean(pon(:,(k-1)*bs+1:k*bs)',1)';
                            d22on(k,:)=mean(d2on((k-1)*bs+1:k*bs,:),1);
                        end
                        for k=1:floor(size(poff,2)/bs)
                            poff2(:,k)=mean(poff(:,(k-1)*bs+1:k*bs)',1)';
                            d22off(k,:)=mean(d2off((k-1)*bs+1:k*bs,:),1);     % average of bs points (bs points)
                        end
                        
                        d22off(abs(d22off)<1e-3) = 0;
                        d22on(abs(d22on)<1e-3) = 0;
                        poff2(abs(poff2)<1e-3) = 0;
                        pon2(abs(pon2)<1e-3) = 0;
                        
                        %Apply PCA to HLS
                        [PCA2,~,PCHLSComp{i}{m}{flyN}{l},reconstructed] = PCA_decode([poff2';pon2']);
                        poff2PCA = PCA2(1:end/2,:);
                        pon2PCA = PCA2(end/2+1:end,:);
                        
                        HLSDistRecOn{flyN,l} = mean(reconstructed(1:end/2,:),1);
                        HLSDistRecOn{flyN,l} = HLSDistRecOn{flyN,l}./sum(HLSDistRecOn{flyN,l});
                        HLSDistRecOff{flyN,l} = mean(reconstructed(end/2+1:end,:),1);
                        HLSDistRecOff{flyN,l} = HLSDistRecOff{flyN,l}./sum(HLSDistRecOff{flyN,l});
                        
                        % raw per fly
                        type{i}.fly{flyN}.rawX{l}=[;d22off;d22on];                          % raw speed and curvature averaged over m points
                        type{i}.fly{flyN}.rawY{l}=[zeros(size(d22off,1),1);ones(size(d22on,1),1)];  % 0/1 for on off speed and curvature
                        % HLS per fly
                        type{i}.fly{flyN}.X{l}=[poff2PCA;pon2PCA];
                        type{i}.fly{flyN}.Y{l}=[zeros(size(poff2PCA,1),1);ones(size(pon2PCA,1),1)];
                        
                        % raw for all flies
                        rawX{l}{i}=[rawX{l}{i};d22off;d22on];
                        rawY{l}{i}=[rawY{l}{i};zeros(size(d22off,1),1);ones(size(d22on,1),1)];
                        % HLS for all flies
                        BX{l}{i}=[BX{l}{i};poff2';pon2'];
                        BY{l}{i}=[BY{l}{i};zeros(size(poff2',1),1);ones(size(pon2',1),1)];
                        
                        display([num2str(flyN) ':' num2str(l)])
                        
                        
                        % fit a logit regression model for observed responses rawY on the predictors rawX (speed/curvature)
                        [type{i}.fly{flyN}.rawB{l},~,~] = glmfit(type{i}.fly{flyN}.rawX{l},type{i}.fly{flyN}.rawY{l},'binomial','link','logit','Options',opts);
                        warnArray.ind{flyN}.raw{l} = checkwarning();
                        % Compute the predicted values for the generalized linear model
                        type{i}.fly{flyN}.rawYhat{l}=glmval(type{i}.fly{flyN}.rawB{l},type{i}.fly{flyN}.rawX{l},'logit');
                        % Compute the percentage of correct predicted values (max predictive power)
                        type{i}.fly{flyN}.PCraw{m}{l} = max(sum((type{i}.fly{flyN}.rawYhat{l}>dc)==repmat(type{i}.fly{flyN}.rawY{l},1,length(dc)))/length(type{i}.fly{flyN}.rawY{l}));
                        
                        % fit a logit regression model for High Level States
                        [type{i}.fly{flyN}.B{l}.p{m},dev,stats] = glmfit(type{i}.fly{flyN}.X{l},type{i}.fly{flyN}.Y{l},'binomial','link','logit','Options',opts);
                        warnArray.ind{flyN}.HLS{l} = checkwarning();
                        type{i}.fly{flyN}.Yhat{l}=glmval(type{i}.fly{flyN}.B{l}.p{m},type{i}.fly{flyN}.X{l},'logit');
                        type{i}.fly{flyN}.PCHLS{m}{l} = max(sum((type{i}.fly{flyN}.Yhat{l}>dc)==repmat(type{i}.fly{flyN}.Y{l},1,length(dc)))/length(type{i}.fly{flyN}.Y{l}));
                        
                    end
                end
            end
            
            for l = 1:2
                [BX{l}{i},~,BXHLSComp{i}{m}{l},~] = PCA_decode(BX{l}{i});
                
                [Braw{l}{i},dev,stats] = glmfit(rawX{l}{i},rawY{l}{i},'binomial','link','logit','Options',opts);
                rawYhat{l}{i}=glmval(Braw{l}{i},rawX{l}{i},'logit');
                [rawpercentCorrect{l}(i),loc]=max(sum(((rawYhat{l}{i}>=dc) == repmat(rawY{l}{i},1,length(dc))))/length(rawY{l}{i}));
                warnArray.group{flyN}.raw{35} = checkwarning();
                
                [BHLS{m}{l}{i},dev,stats] = glmfit(BX{l}{i},BY{l}{i},'binomial','link','logit','Options',opts);
                Yhat{l}{i}=glmval(BHLS{m}{l}{i},BX{l}{i},'logit');
                [HLSpercentCorrect{l}(i),loc]=max(sum(((Yhat{l}{i}>=dc) == repmat(BY{l}{i},1,length(dc))))/length(BY{l}{i}));
                warnArray.group{flyN}.HLS{35} = checkwarning();
                
                %Apply predictor to each fly individually
                for flyN=idx2
                    if sum(flyN == bad_Fly_List)<1
                        type{i}.fly{flyN}.rawYhatALL{m}{l} = glmval(Braw{l}{i},type{i}.fly{flyN}.rawX{l},'logit');
                        type{i}.fly{flyN}.PCrawALL{m}{l} = sum((type{i}.fly{flyN}.rawYhatALL{m}{l}>1/2)==type{i}.fly{flyN}.rawY{l})/length(type{i}.fly{flyN}.rawY{l});
                        warnArray.group{flyN}.raw{flyN} = checkwarning();
                        
                        if size(type{i}.fly{flyN}.X{l},2)>=size(BHLS{m}{l}{i},1)
                            type{i}.fly{flyN}.X{l} = type{i}.fly{flyN}.X{l}(:,1:size(BHLS{m}{l}{i},1)-1);
                        end
                        
                        BHLS_fly = BHLS{m}{l}{i}(1:size(type{i}.fly{flyN}.X{l},2)+1);
                        type{i}.fly{flyN}.HLSYhatALL{m}{l} = glmval(BHLS_fly,type{i}.fly{flyN}.X{l},'logit');
                        type{i}.fly{flyN}.PCHLSALL{m}{l} = sum((type{i}.fly{flyN}.HLSYhatALL{m}{l}>1/2)==type{i}.fly{flyN}.Y{l})/length(type{i}.fly{flyN}.Y{l});
                        warnArray.group{flyN}.HLS{flyN} = checkwarning();
                        
                    end
                end
            end
        end
    end
end

end



