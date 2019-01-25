function [meanDist,meanDist2,sampDist,sampDist2,Z1,aStruct,Dn,otherInfo] = ...
    SyntheticFlyDecoder(data,model,likely_high_state_sorted,NumPoints,...
    durVec,B,NumPointsAll,flyClustNdx,options,saveDat)
%SyntheticFlyDecoder Takes in information about the raw distribution of HLS
%used in the GLM analysis, the duration the each fly spends
%inside/outside the odor zone,before/after application of the odor in order
%to generate synthetic tracks and conduct analysis on these tracks

% Inputs:
%    NumPoints: 1x2 cell {1} = During, {2} = Before, each cell is a nx2
%    matrix where (:,1) = inside, (:,2) = outside, n = # of flies
%    durVec: 2x2 cell where rows are for location (in,out), columns are for
%    state (During,Before). Each cell is 1x10 cell, one for each HLS. Each
%    HLS cell is composed of 1xn double where n is max cont. len of HLS
%    B: 2x3 cell where rows are location (in,out), columns are for state 
%    (During,Before,Difference). Each cell contains n HLS distributions
%    NumPointsAll: 1x2 cell {1} = During, {2} = Before, each cell is a 34x2
%    matrix where (:,1) = inside, (:,2) = outside
%    saveDat: 1 = yes, 0 = no

% Outputs:
%    meanDist: Mean weighted HLS distribution for all emp flies used to
%    generate sampDist
%    meanDist2: Mean weighted HLS track length distribution for all emp 
%    flies used to generate sampDist2
%    sampDist: Distribution used as first order library to generate
%    synthetic tracks
%    sampDist2: Distribution used as second order library to generate
%    synthetic tracks
%    Z1: Contains dendrogram information for one synthetic fly run (34 flies)
%    aStruct: Matrix containing information about HLS distributions
%    Dn: Matrix containing information about euclidean and Djs distances
%    otherInfo: misc info used in the generating the synthetic flies

% Note that all inputs are generated from GLM_Decode. Before and during
% indexes are switched in this analysis as we designate before then during
%
% 2017, Liangyu Tao

%close all
% switch before and during so that it goes before,during,diff
B = B(~cellfun(@isempty,B));
NumPoints_new = NumPoints([2,1]);                                           % same as above
NumPoints_All = NumPointsAll([2,1]);
durVec_new = durVec([2,1;4,3]);
meanDist = cell(1,length(flyClustNdx));
meanDist2 = cell(1,length(flyClustNdx));
sampDist = cell(1,length(flyClustNdx));
sampDist2 = cell(1,length(flyClustNdx));
Z1 = cell(1,length(flyClustNdx));
aStruct = cell(1,length(flyClustNdx));
Dn = cell(1,length(flyClustNdx));
flyClustNdx = flyClustNdx(~cellfun(@isempty,flyClustNdx));
NumTrials = options.numSynth;                                                            % number of simulations

for c = 1:length(flyClustNdx)
    if ~isempty(flyClustNdx{c})
        B_new{c} = B{c}(:,[2 1 3]);                                                       % before: B = {:,1}, during B = {:,2}
        flies = flyClustNdx{c};
        
        leng = [];
        medLeng = zeros(2,2);
        for state = 1:2                                             % 1 = before, 2 = during
            for loc = 1:2                                           % 1 = inside, 2 = outside
                medLeng(state,loc) = round(median(NumPoints_All{state}(:,loc)));
                leng = [leng medLeng(state,loc)];
            end
        end
        otherInfo{c} = [];
        otherInfo{c}.medLeng = medLeng;
        otherInfo{c}.numPoints = NumPoints_All;
        otherInfo{c}.durVec = durVec_new;
        otherInfo{c}.B_new = B_new{c};

        leng_to_consider = 1:4;
        numFlies = length(flies);
        
        [meanDist{c},meanDist2{c},sampDist{c},sampDist2{c}] = createLibrary(B_new{c},NumPoints_new,durVec_new,flies);
%        [Z1{c},aStruct{c}] = generateSyntheticTrack(sampDist{c},sampDist2{c},NumTrials,leng,numFlies);
        
        [TPAll,pi0] = generateTP(data,model,likely_high_state_sorted,flies);
        [Z1{c},aStruct{c}] = generateSyntheticTrackTP(TPAll,pi0,NumTrials,leng,numFlies);

        Dn{c} = calcDistances(aStruct{c},B_new{c},NumTrials,leng,numFlies);
    end
end

if saveDat == 1
    save(options.synthFile,'-v7.3')
end

end

function [meanDist,meanDist2,sampDist,sampDist2] = createLibrary(B,NumPoints,durVec,flies)
% close all
HLNum = size(durVec{1},2);
meanDist2 = cell(2,2,HLNum);
sampDist2 = cell(2,2,HLNum);                                                   % library for second order based on durations
for state = 1:2
    for loc = 1:2
        for HLS = 1:HLNum
            tempLength = sort(durVec{state,loc}{HLS});
            if ~isempty(tempLength)
                meanDist2{state,loc,HLS} = round(1000*histcounts(tempLength,[0:max(tempLength)+1])/length(tempLength));
                for i = 1:length(meanDist2{state,loc,HLS})
                    sampDist2{state,loc,HLS} = [sampDist2{state,loc,HLS}, i*ones(1,meanDist2{state,loc,HLS}(i))];
                end
                sampDist2{state,loc,HLS}(end+1:1000)=0;
                sampDist2{state,loc,HLS}=sampDist2{state,loc,HLS}(end-999:end);
                sampDist2{state,loc,HLS} = sort(sampDist2{state,loc,HLS});
            end
        end
    end
end
% generate mean distribution and create a library based on this avg dist
meanDist = cell(2,2);
sampDist = cell(2,2);                                                       % library for first order based on probability
for state = 1:2                                             % 1 = before, 2 = during
    for loc = 1:2                                           % 1 = inside, 2 = outside
        Freq = repmat(NumPoints{state}(flies,loc),1,HLNum)';
        meanDist{state,loc} = round(sum(B{loc,state}.*Freq,2)./sum(NumPoints{state}(flies,loc)),3);
        for HLS = 1:HLNum
            sampDist{state,loc} = [sampDist{state,loc}, HLS*ones(1,1000*meanDist{state,loc}(HLS))];
        end
        sampDist{state,loc}(end+1:1000)=1;
        if sum(meanDist{state,loc}) ~= 1
            meanDist{state,loc} = meanDist{state,loc}./sum(meanDist{state,loc});
        end
        sum(meanDist{state,loc})
    end
end
end

function [TPAll,pi0] = generateTP(data,model,likely_high_state_sorted,flies)
HLNum = model.Qdim;
len = length(data{1}(1,:));
TPAll = cell(2);%pi0 = cell(2);
for i = 1:4
    TPAll{i} = zeros(HLNum,HLNum);
end
for i = flies
    TP = zeros(model.Qdim,model.Qdim,len-1);
    for t = 1:len-1
        if likely_high_state_sorted(i,t)>0 && likely_high_state_sorted(i,t+1)>0
            TP(likely_high_state_sorted(i,t),likely_high_state_sorted(i,t+1),t) ...
                = 1+TP(likely_high_state_sorted(i,t),likely_high_state_sorted(i,t+1),t);
        end
    end
    for state = 1:2
        for loc = 1:2
            c = data{i}(11,1:end-1)==(state-1)&(1-data{i}(12,1:end-1)==(loc-1));
            TPAll{state,loc} = TPAll{state,loc}+sum(TP(:,:,c),3);
        end
    end
end

pi0 = cell(2,2);
for state = 1:2
    for loc = 1:2
        TP_sorted=TPAll{state,loc};
        TPTmp = TP_sorted./repmat(sum(TP_sorted,2),1,model.Qdim);
        TPAll{state,loc} = TPTmp;
        pi0{state,loc} = sum(TP_sorted,2)'./sum(TP_sorted(:));
        %pi0{state,loc} = ones(1,size(TP_sorted,1))./size(TP_sorted,1);
    end
end

end

% based on duration distributions
function [Z1,aStruct] = generateSyntheticTrack(sampDist,sampDist2,nIt,leng,numFlies)
HLNum = size(sampDist2,3);
dist = cell(2,2);
meanDistNew = cell(2,2,4);
aStruct = [];
for numTrials = 1:nIt
    for i = 1:length(leng)
        % create a synthetic HLS track based on the library
        synth_HLS_All_Flies = cell(2,2);
        %len = 140;
        len = leng(i);
        for state = 1:2                                             % 1 = before, 2 = during
            for loc = 1:2                                           % 1 = inside, 2 = outside
                synth = round(1000*rand(numFlies,len));
                synth(synth==0) = 1;
                synth_HLS_All_Flies{state,loc} = sampDist{state,loc}(synth);
            end
        end
        
        synTrackLength = cell(2,2);
        synthTrackAll = cell(2,2);
        for state = 1:2                                             % 1 = before, 2 = during
            for loc = 1:2                                           % 1 = inside, 2 = outside
                tempTrackLength = synth_HLS_All_Flies{state,loc};
                synthTrackAll{state,loc} = zeros(numFlies,len);
                for j = 1:len
                    for HLS = 1:HLNum
                        ndx = find(synth_HLS_All_Flies{state,loc}(:,j)==HLS);
                        syn = round(1000*rand(length(ndx),1));
                        syn(syn==0) = 1;
                        tempTrackLength(ndx,j) = sampDist2{state,loc,HLS}(syn);
                    end
                end
                
                synTrackLength{state,loc} = tempTrackLength;
                
                for fly = 1:numFlies
                    initState = synth_HLS_All_Flies{state,loc}(fly,:);
                    StateLength = tempTrackLength(fly,:);
                    synthTrackTemp = [];
                    idx = 1;
                    while length(synthTrackTemp) <=len
                        synthTrackTemp = [synthTrackTemp, initState(idx)*ones(1,StateLength(idx)+1)];
                        idx = idx+1;
                    end
                    synthTrackAll{state,loc}(fly,:) = synthTrackTemp(1:len);
                end
            end
        end
        aStruct.synthTrackAll{i,numTrials} = synthTrackAll;
        
        % get the distribution based on the synthetic HLS track
        synth_HLS_Dist = cell(2,2);
        for state = 1:2                                             % 1 = before, 2 = during
            for loc = 1:2                                           % 1 = inside, 2 = outside
                for fly = 1:numFlies
                    for hls = 1:HLNum
                        synth_HLS_Dist{state,loc}(fly,hls) = sum(synthTrackAll{state,loc}(fly,:)==hls);
                    end
                    synth_HLS_Dist{state,loc}(fly,:) = synth_HLS_Dist{state,loc}(fly,:)./sum(synth_HLS_Dist{state,loc}(fly,:));
                end
                meanDistNew{state,loc,i}(numTrials,:) = mean(synth_HLS_Dist{state,loc});
                aStruct.SynthDist{state,loc}{i,numTrials} = synth_HLS_Dist{state,loc};
            end
        end
        
        fly = 1:1:numFlies;
        % plot dendrogram for these flies
        Z1 = cell(2,2);
        s = {'before', 'during'};
        l = {'inside', 'outside'};
%         if numTrials == 1
%             figure
%         end
%         kk = 1;
        for loc = 1:2
            for state = 1:2
                Z1{state,loc} = linkage(synth_HLS_Dist{loc,state},'ward');
                flylabel = cellfun(@num2str, num2cell(fly), 'UniformOutput', false);
                dist{state,loc}(i,numTrials) = Z1{state,loc}(end);
%                 if numTrials == 1
%                     subplot(2,2,kk);dendrogram(Z1{state,loc},0,'labels',flylabel);
%                     ylim([0 ceil(100*Z1{1,1}(end))/50])
%                     title([s{state} ' ' l{loc}])
%                     kk = kk+1;
%                 end
                
            end
        end
%         if numTrials == 1
%             suptitle([num2str(len) ' samples'])
%             set(gcf,'position',[9 49 800 918])
%             %print('-dpsc2',[['PDF_Files/Dendrograms_Synthetic'] '.ps'],'-loose','-append');
%         end
    end
    if mod(numTrials,10)==0
        display(numTrials);
    end
end
aStruct.meanDistNew = meanDistNew;
aStruct.maxDistance = dist;
end

% Based on transition probabilities
function [Z1,aStruct] = generateSyntheticTrackTP(TPAll,pi0,nIt,leng,numFlies)
HLNum = length(pi0{1});
dist = cell(2,2);
meanDistNew = cell(2,2,4);
aStruct = [];
for numTrials = 1:nIt
    for i = 1:length(leng)
        len = leng(i);
        synthTrackAll = cell(2);
        % create a synthetic HLS track based on the library
        for state = 1:2
            for loc = 1:2
                initState = discretesample(pi0{state,loc},numFlies);
                [synthTrackAll{state,loc}] = syntheticFlyGeneration(initState,TPAll{state,loc},len,HLNum);
            end
        end
        aStruct.synthTrackAll{i,numTrials} = synthTrackAll;
        
        % get the distribution based on the synthetic HLS track
        synth_HLS_Dist = cell(2,2);
        for state = 1:2                                             % 1 = before, 2 = during
            for loc = 1:2                                           % 1 = inside, 2 = outside
                for fly = 1:numFlies
                    for hls = 1:HLNum
                        synth_HLS_Dist{state,loc}(fly,hls) = sum(synthTrackAll{state,loc}(fly,:)==hls);
                    end
                    synth_HLS_Dist{state,loc}(fly,:) = synth_HLS_Dist{state,loc}(fly,:)./sum(synth_HLS_Dist{state,loc}(fly,:));
                end
                meanDistNew{state,loc,i}(numTrials,:) = mean(synth_HLS_Dist{state,loc});
                aStruct.SynthDist{state,loc}{i,numTrials} = synth_HLS_Dist{state,loc};
            end
        end
        
        % dendrogram analysis (optional)
        fly = 1:1:numFlies;
        Z1 = cell(2,2);
        for loc = 1:2
            for state = 1:2
                Z1{state,loc} = linkage(synth_HLS_Dist{state,loc},'ward');
%                 flylabel = cellfun(@num2str, num2cell(fly), 'UniformOutput', false);
                dist{state,loc}(i,numTrials) = Z1{state,loc}(end);
            end
        end
    end
    if mod(numTrials,10)==0
        display(numTrials);
    end
end
aStruct.meanDistNew = meanDistNew;
aStruct.maxDistance = dist;

end

function [HLSSyn] = syntheticFlyGeneration(initState,TP,len,HLNum)
%generate HLS distributions
HLSDist = cell(1,HLNum);
TP(isnan(TP)) = 1/HLNum;
for i = 1:HLNum
    for j = 1:HLNum
        HLSDist{i} = [HLSDist{i} j.*ones(1,round(TP(i,j).*10000))];
    end
    HLSDist{i}(end+1:10000)=2;
end

%generate tracks of HLS
HLSSyn = zeros(length(initState),len);HLSSyn(:,1) = initState;
for i = 1:len
    for HLS = 1:HLNum
        currHLS = (HLSSyn(:,i)==HLS);
        nextHLSNdx = ceil(10000.*rand(sum(currHLS),1));
        HLSSyn(currHLS,i+1) = HLSDist{HLS}(nextHLSNdx);
    end
end
end

