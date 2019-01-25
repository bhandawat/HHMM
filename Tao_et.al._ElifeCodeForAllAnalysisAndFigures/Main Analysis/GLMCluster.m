function [NumPoints,durVec,PCrawALL,PCHLSALL,All,Ind] = GLMCluster(data,model,likely_high_state_sorted,kMeans,c2cons)
%GLM_decode Takes in the HLS for each fly track to fit to a logistic
%regression model as a means of decoding Fly behavior. Dendrogram analysis
%is further conducted here

% Inputs:
%    likely_high_state_by_fly: 34x10799 matrix of unsorted HLS tracks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) Transpose Property Matrices
%       2.) Speed and curvature for HLS
%       3.) color scheme for each cluster
%       4.) Speed and Curvature for LL states
%       5.) All tracks with initial velocity going upwards
%       6.) Decoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that the output of this will go into other analysis
%
% 2017, Bhandawat Lab

%clearvars -except kMeans
%clc
dc=0:0.01:1;
opts=statset('glmfit');
opts.MaxIter=1000;
% number of points to average over
bs = 30;
%warning('off', 'stats:pca:ColRankDefX')
warning('off')

lastwarn('')
warnArray = [];
warnArray.group = cell(1,35);
warnArray.ind = cell(1,34);
bad_points = [];
bad_points.j = [];
bad_points.l = [];
NumPoints = cell(1,2);
data = post_hoc_preprocess_LT(data);

durVec = cell(2,2);
for l = 1:4
    durVec{l}=cell(1,model.Qdim);
end

flysToConsider = cell(2,10);
Ind = [];
All = [];
rawpercentCorrect = [];
HLSpercentCorrect = [];
PCrawALL = [];
PCHLSALL = [];
bad_Fly_List = [];

for i=1:length(c2cons)                                                           % number of clusters based on HHMM model
    idx2=c2cons{i};  % idx of datapoints assigned to cluster i
    np=length(idx2);
    
    if(np>0)
        for locK = 1:2                                                          % 1 = inside, 2 = outside
            numCluster = kMeans.data{locK,3}.optCluster;
            
            for clusters = 1:numCluster
                flysToConsider{locK,clusters} = kMeans.data{locK,3}.cluster{numCluster}{clusters};
            end
            
            for clusters = 1:numCluster
                f = 0;
                for flys = flysToConsider{locK,clusters}
                    f = f+1;
                    d1=data{idx2(flys)};
                    p2 = zeros(model.Qdim,size(d1,2));
                    idxon = cell(1,2);idxoff = cell(1,2);
                    for k=1:model.Qdim
                        
                        ndx = find(likely_high_state_sorted(idx2(flys),:)==k);
                        p2(k,ndx) = 1;
                        idxoff{2}=find(d1(11,1:end-1)<1 & d1(12,1:end-1)<1);    % odor off and fly outside
                        idxoff{1}=find(d1(11,1:end-1)<1 & d1(12,1:end-1)==1);   % odor off and fly inside
                        idxon{2}=find(d1(11,1:end-1)==1 & d1(12,1:end-1)<1);    % odor on and fly outside
                        idxon{1}=find(d1(11,1:end-1)==1 & d1(12,1:end-1)==1);   % odor on and fly inside
                    end
                    
                    if length(idxon{1})<91 || length(idxoff{1})<91 || length(idxon{2})<91 || length(idxoff{2})<91
                        display(['Fly ' num2str(idx2(flys)) ' is a bad fly'])
                        display([num2str(length(idxon{1})) '  ' num2str(length(idxoff{1}))...
                            '  ' num2str(length(idxon{2})) '  ' num2str(length(idxoff{2}))])
                        bad_Fly_List = [bad_Fly_List, idx2(flys)];
                    else
                        for loc = 1:2                                               % 1 = inside, 2 = outside
                            
                            %setting chance to zero by tossing out points until
                            %idxon and idxoff are same length
                            pon=p2(:,idxon{loc});                                 % Prob of HL states given stim on
                            poff=p2(:,idxoff{loc});                               % Prob of HL states given stim off
                            ndx2 = sum(pon,1)==0;
                            idxon{loc}(ndx2) = [];
                            pon(:,ndx2) = [];
                            ndx2 = sum(poff,1)==0;
                            idxoff{loc}(ndx2) = [];
                            poff(:,ndx2) = [];
                            len=min(size(pon,2),size(poff,2));
                            pon=pon(:,1:len);
                            poff=poff(:,1:len);
                            
                            idxoff{loc}=idxoff{loc}(1:len);
                            idxon{loc}=idxon{loc}(1:len);
                            
                            d2on=d1(7:8,idxon{loc})';
                            d2on(8,:)=abs(d2on(8,:));                           % speed and curvature given stim on
                            d2off=d1(7:8,idxoff{loc})';
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
                            [PCA2,~,~,reconstructed] = PCA_decode([poff2';pon2']);
                            poff2PCA = PCA2(1:end/2,:);
                            pon2PCA = PCA2(end/2+1:end,:);
                            
                            % raw per fly
                            Ind{i}.KC{locK,clusters}.fly{f}.rawX{loc} = [d22off;d22on];
                            Ind{i}.KC{locK,clusters}.fly{f}.rawY{loc} = [zeros(size(d22off,1),1);ones(size(d22on,1),1)];
                            tempRawX = Ind{i}.KC{locK,clusters}.fly{f}.rawX{loc};
                            tempRawY = Ind{i}.KC{locK,clusters}.fly{f}.rawY{loc};
                            % HLS per fly
                            Ind{i}.KC{locK,clusters}.fly{f}.X{loc}=[poff2PCA;pon2PCA];
                            Ind{i}.KC{locK,clusters}.fly{f}.Y{loc}=[zeros(size(poff2PCA,1),1);ones(size(pon2PCA,1),1)];
                            tempHLSX = Ind{i}.KC{locK,clusters}.fly{f}.X{loc};
                            tempHLSY = Ind{i}.KC{locK,clusters}.fly{f}.Y{loc};
                            
                            % raw for all flies
                            if f == 1
                                All{i}.KC{locK,clusters}.rawX{loc}=[d22off;d22on];
                                All{i}.KC{locK,clusters}.rawY{loc}=[zeros(size(d22off,1),1);ones(size(d22on,1),1)];
                            else
                                All{i}.KC{locK,clusters}.rawX{loc}=[All{i}.KC{locK,clusters}.rawX{loc};d22off;d22on];
                                All{i}.KC{locK,clusters}.rawY{loc}=[All{i}.KC{locK,clusters}.rawY{loc};zeros(size(d22off,1),1);ones(size(d22on,1),1)];
                            end
                            % HLS for all flies
                            if f == 1
                                All{i}.KC{locK,clusters}.X{loc}=[poff2';pon2'];
                                All{i}.KC{locK,clusters}.Y{loc}=[zeros(size(poff2',1),1);ones(size(pon2',1),1)];
                                All{i}.KC{locK,clusters}.Cutoff{loc} = [1, 1+size(poff2,2)*2];
                            else
                                All{i}.KC{locK,clusters}.X{loc}=[All{i}.KC{locK,clusters}.X{loc};poff2';pon2'];
                                All{i}.KC{locK,clusters}.Y{loc}=[All{i}.KC{locK,clusters}.Y{loc};zeros(size(poff2',1),1);ones(size(pon2',1),1)];
                                All{i}.KC{locK,clusters}.Cutoff{loc} = [All{i}.KC{locK,clusters}.Cutoff{loc}, All{i}.KC{locK,clusters}.Cutoff{loc}(end)+size(poff2,2)*2];
                            end
                            
                            % fit a logit regression model for observed responses rawY on the predictors rawX (speed/curvature)
                            [Ind{i}.KC{locK,clusters}.fly{f}.rawB{loc},~,~] = glmfit(tempRawX,tempRawY,'binomial','link','logit','Options',opts);
                            Ind{i}.KC{locK,clusters}.fly{f}.rawYhat{loc}=glmval(Ind{i}.KC{locK,clusters}.fly{f}.rawB{loc},tempRawX,'logit');
                            Ind{i}.KC{locK,clusters}.fly{f}.PCraw{loc} = max(sum((Ind{i}.KC{locK,clusters}.fly{f}.rawYhat{loc}>dc)==repmat(tempRawY,1,length(dc)))/length(tempRawY));
                            
                            % fit a logit regression model for High Level States
                            [Ind{i}.KC{locK,clusters}.fly{f}.B{loc},~,~] = glmfit(tempHLSX,tempHLSY,'binomial','link','logit','Options',opts);
                            Ind{i}.KC{locK,clusters}.fly{f}.Yhat{loc}=glmval(Ind{i}.KC{locK,clusters}.fly{f}.B{loc},tempHLSX,'logit');
                            Ind{i}.KC{locK,clusters}.fly{f}.PC{loc} = max(sum((Ind{i}.KC{locK,clusters}.fly{f}.Yhat{loc}>dc)==repmat(tempHLSY,1,length(dc)))/length(tempHLSY));
                        end
                    end
                end
                
                for loc = 1:2
                    rawX = All{i}.KC{locK,clusters}.rawX{loc};
                    rawY = All{i}.KC{locK,clusters}.rawY{loc};
                    [HLSX,~,~,~] = PCA_decode(All{i}.KC{locK,clusters}.X{loc});
                    HLSY = All{i}.KC{locK,clusters}.Y{loc};
                    
                    All{i}.KC{locK,clusters}.PCAX{loc} = HLSX;
                    
                    [Braw,dev,stats] = glmfit(rawX,rawY,'binomial','link','logit','Options',opts);
                    rawYhat=glmval(Braw,rawX,'logit');
                    [rawpercentCorrect{i,loc}.KC{locK,clusters},loc2]=max(sum(((rawYhat>=dc) == repmat(rawY,1,length(dc))))/length(rawY));
                    %warnArray.group{j}.raw{35} = checkwarning();
                    
                    [BHLS,dev,stats] = glmfit(HLSX,HLSY,'binomial','link','logit','Options',opts);
                    Yhat=glmval(BHLS,HLSX,'logit');
                    [HLSpercentCorrect{i,loc}.KC{locK,clusters},loc2]=max(sum(((Yhat>=dc) == repmat(HLSY,1,length(dc))))/length(HLSY));
                    %warnArray.group{j}.HLS{35} = checkwarning();
                    
                    %Apply predictor to each fly individually
                    for j=1:f
                        if sum(j == bad_Fly_List)<1
                            ndxTemp = All{i}.KC{locK,clusters}.Cutoff{loc}(j):All{i}.KC{locK,clusters}.Cutoff{loc}(j+1)-1;
                            RawFlyX = rawX(ndxTemp,:);
                            RawFlyY = rawY(ndxTemp,:);
                            HLSFlyX = HLSX(ndxTemp,:);
                            HLSFlyY = HLSY(ndxTemp,:);
                            
                            rawYhatALL=glmval(Braw,RawFlyX,'logit');
                            PCrawALL{i,loc}.KC{locK,clusters}(j) = sum((rawYhatALL>1/2)==RawFlyY)/length(RawFlyY);
                            HLSYhatALL=glmval(BHLS,HLSFlyX,'logit');
                            PCHLSALL{i,loc}.KC{locK,clusters}(j) = sum((HLSYhatALL>1/2)==HLSFlyY)/length(HLSFlyY);
                        end
                    end
                end
                
            end
            
            for loc = 1:2
                rawX = All{i}.KC{locK,1}.rawX{loc};
                rawY = All{i}.KC{locK,1}.rawY{loc};
                HLSX = All{i}.KC{locK,1}.X{loc};
                HLSY = All{i}.KC{locK,1}.Y{loc};
                cNdx = zeros(1,7);
                cNdx(2) = length(All{i}.KC{locK,1}.X{loc});
                for clusters = 1:numCluster
                    rawX = [rawX;All{i}.KC{locK,clusters}.rawX{loc}];
                    rawY = [rawY;All{i}.KC{locK,clusters}.rawY{loc}];
                    HLSX = [HLSX;All{i}.KC{locK,clusters}.X{loc}];
                    HLSY = [HLSY;All{i}.KC{locK,clusters}.Y{loc}];
                    cNdx(clusters+1) = cNdx(clusters)+length(All{i}.KC{locK,clusters}.X{loc});
                end
                
                [HLSX,~,~,~] = PCA_decode(HLSX);
                All{i}.KC{locK,clusters}.PCAX{loc} = HLSX;
                
                [Braw,dev,stats] = glmfit(rawX,rawY,'binomial','link','logit','Options',opts);
                [BHLS,dev,stats] = glmfit(HLSX,HLSY,'binomial','link','logit','Options',opts);
                for clusters = 1:numCluster
                    
                    %Apply predictor to each fly individually
                    for j=1:length(flysToConsider{locK,clusters})
                        if sum(j == bad_Fly_List)<1
                            ndxTemp = All{i}.KC{locK,clusters}.Cutoff{loc}(j)+cNdx(clusters):All{i}.KC{locK,clusters}.Cutoff{loc}(j+1)-1+cNdx(clusters);
                            RawFlyX = rawX(ndxTemp,:);
                            RawFlyY = rawY(ndxTemp,:);
                            HLSFlyX = HLSX(ndxTemp,:);
                            HLSFlyY = HLSY(ndxTemp,:);

                            rawYhatALL=glmval(Braw,RawFlyX,'logit');
                            PCrawALL{i,loc}.AC{locK,clusters}(j) = sum((rawYhatALL>1/2)==RawFlyY)/length(RawFlyY);
                            HLSYhatALL=glmval(BHLS,HLSFlyX,'logit');
                            PCHLSALL{i,loc}.AC{locK,clusters}(j) = sum((HLSYhatALL>1/2)==HLSFlyY)/length(HLSFlyY);
                        end
                    end
                end
            end
            
            
        end
    end
end

end






