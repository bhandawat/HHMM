function [Dn] = calcDistances(aStruct,B,numTrials,leng,numFlies)
%calcDistances takes in information about the synthetic fly tracks and
%computes the euclidean distances and Jensen Shannon Divergences between
%flies.

% Inputs:
%    aStruct: Matrix containing information about HLS distributions
%    B: 2x3 cell where rows are location (in,out), columns are for state 
%    (Before,During,Difference). Each cell contains n HLS distributions
%    numTrials: Number of iterations of synth tracks generated
%    leng: list of median length (state-1)*2+loc
%    numFlies: Number of flies to analysize (cluster size)

% Outputs:
%    Dn: Matrix containing information about euclidean and Djs distances
%    otherInfo: misc info used in the generating the synthetic flies

%
% 2017, Liangyu Tao

scenario = {'b_i','b_o';'d_i','d_o'};
HLNum = size(B{1},1);
Dn.scenario = scenario;
% calculate shannon divergence and euclidean distance for emperical
uniqueCompJS = cell(size(B,2)-1,size(B,1));
uniqueCompIJ = cell(size(B,2)-1,size(B,1));
Meanij = cell(size(B,2)-1,size(B,1));
Meanjs = cell(size(B,2)-1,size(B,1));
Medianij = cell(size(B,2)-1,size(B,1));
Medianjs = cell(size(B,2)-1,size(B,1));
NjsTemp = cell(size(B,2)-1,size(B,1));
NijTemp = cell(size(B,2)-1,size(B,1));
DijTemp = cell(size(B,2)-1,size(B,1));
for state = 1:2                                                             % B(state = 2) = before,B(state = 1) = during
    Djs = cell(1,size(B,1));
    Dij = cell(1,size(B,1));
    for loc = 1:size(B,1)
        distTemp = B{loc,state}';
        % add in a small delta to locations where P = 0;
        for j = 1:numFlies
            distTemp2 = distTemp(j,:);
            ndx = distTemp2==0;
            distTemp2(ndx) = distTemp2(ndx)+0.0001;                         % epsilon = 0.0001
            distTemp2(~ndx) = distTemp2(~ndx)-0.0001*sum(ndx)/(HLNum-sum(ndx));
            distTemp(j,:) = distTemp2;
        end
        for j = 1:numFlies
            for k = 1:numFlies
                % calculate shannon divergence
                Djs{loc}(j,k) = jensen_shannon_divergence(distTemp(j,:),distTemp(k,:));
                % calculate euclidean distance
                Dij{loc}(j,k) = sqrt(sum((distTemp(j,:)-distTemp(k,:)).^2));
            end
        end
        DijTemp{state,loc} = Dij{loc};
        
        uniqueCompJS{state,loc} = squareform(Djs{loc});
        Meanjs{state,loc} = mean(uniqueCompJS{state,loc});
        Medianjs{state,loc} = median(uniqueCompJS{state,loc});
        NjsTemp{state,loc} = histcounts(uniqueCompJS{state,loc},[0:0.05:1])/length(uniqueCompJS{state,loc});
        
        uniqueCompIJ{state,loc} = squareform(Dij{loc});
        Meanij{state,loc} = mean(uniqueCompIJ{state,loc});
        Medianij{state,loc} = median(uniqueCompIJ{state,loc});
        NijTemp{state,loc} = histcounts(uniqueCompIJ{state,loc},[0:0.15:3])/length(uniqueCompIJ{state,loc});
    end
end
Dn.jsEmp.all = NjsTemp;
Dn.jsEmp.mean = Meanjs;
Dn.jsEmp.med = Medianjs;
Dn.ijEmp.all = NijTemp;
Dn.ijEmp.mean = Meanij;
Dn.ijEmp.med = Medianij;
Dn.ijEmp.tot = DijTemp;

% calculate shannon divergence and euclidean distance for synthesis
for state = 1:2                                             % 1 = before, 2 = during
    for loc = 1:2                                           % 1 = inside, 2 = outside
        uniqueCompJS = cell(1,numTrials);
        uniqueCompIJ = cell(1,numTrials);
        NjsTemp = cell(1,numTrials);
        NijTemp = cell(1,numTrials);
        Meanjs = cell(1,numTrials);
        Medianjs = cell(1,numTrials);
        Meanij = cell(1,numTrials);
        Medianij = cell(1,numTrials);
        
        for l = 1:length(leng)
            Djs = cell(1,numTrials);
            Dij = cell(1,numTrials);
            for i = 1:numTrials
                distTemp = aStruct.SynthDist{state,loc}{l,i};
                % add in a small delta to locations where P = 0;
                for j = 1:numFlies
                    distTemp2 = distTemp(j,:);
                    ndx = distTemp2==0;
                    distTemp2(ndx) = distTemp2(ndx)+0.0001;                         % epsilon = 0.0001
                    distTemp2(~ndx) = distTemp2(~ndx)-0.0001*sum(ndx)/(HLNum-sum(ndx));
                    distTemp(j,:) = distTemp2;
                end
                for j = 1:numFlies
                    for k = 1:numFlies
                        % calculate shannon divergence
                        Djs{i}(j,k) = real(jensen_shannon_divergence(distTemp(j,:),distTemp(k,:)));
                        % calculate euclidean distance
                        Dij{i}(j,k) = sqrt(sum((distTemp(j,:)-distTemp(k,:)).^2));
                    end
                end
                uniqueCompJS{l,i} = squareform(Djs{i});
                NjsTemp{l,i} = histcounts(uniqueCompJS{l,i},[0:0.05:1])/length(uniqueCompJS{l,i});
                Meanjs{l,i} = mean(uniqueCompJS{l,i});
                Medianjs{l,i} = median(uniqueCompJS{l,i});
                
                uniqueCompIJ{l,i} = squareform(Dij{i});
                NijTemp{l,i} = histcounts(uniqueCompIJ{l,i},[0:0.15:3])/length(uniqueCompIJ{l,i});
                Meanij{l,i} = mean(uniqueCompIJ{l,i});
                Medianij{l,i} = median(uniqueCompIJ{l,i});
            end
            Dn.ijSyn.(scenario{state,loc}).tot{l} = Dij;
        end
        Dn.jsSyn.(scenario{state,loc}).all = NjsTemp;
        Dn.jsSyn.(scenario{state,loc}).mean = Meanjs;
        Dn.jsSyn.(scenario{state,loc}).med = Medianjs;
        Dn.ijSyn.(scenario{state,loc}).all = NijTemp;
        Dn.ijSyn.(scenario{state,loc}).mean = Meanij;
        Dn.ijSyn.(scenario{state,loc}).med = Medianij;
    end
end
end