function [] = SyntheticClusterDecode(data,model,likely_high_state_sorted,...
    kMeans,B,indDurVec,durVec,NumPoints,c2cons,options)
%SyntheticClusterDecode Takes in information from kMeans on empirical flies
%and synthetic flies based on the median fly to generate synthetic flies
%based on the median of each cluster of flies

% Inputs:
%    kMeans: Cell containing kMeans clustering of flies.
%    B: 2x3 cell where rows are location (in,out), columns are for state 
%    (Before,During,Difference). Each cell contains n HLS distributions
%    indDurVec: N/A
%    durVec: 2x2 cell where rows are for location (in,out), columns are for
%    state (During,Before). Each cell is 1x10 cell, one for each HLS. Each
%    HLS cell is composed of 1xn double where n is max cont. len of HLS
%    NumPoints: 1x2 cell {1} = During, {2} = Before, each cell is a nx2
%    matrix where (:,1) = inside, (:,2) = outside, n = # of flies
%    datName: name of data file to save to

%
% Note that there is no output, but will save results to data file
%
% 2017, Liangyu Tao

S = cell(1,length(c2cons));
l = {'Inside','Outside'};
k = 0;
for cl = 1:length(c2cons)
    if ~isempty(c2cons{cl})
        k = k+1;
        for i = 1:2
            for j = 1:2                                                             %before,during
                optCluster = kMeans{k}.maxCluster.data{i,j}.optCluster;
                for c = 1:optCluster
                    currentFlies = kMeans{k}.maxCluster.data{i,j}.cluster{1,optCluster}{c};
                    B_C = cell(2,3);
                    NumPoints_C = cell(1,2);
                    for loc = 1:2
                        for state = 1:3
                            B_C{loc,state}=B{cl}{loc,state}(:,currentFlies);
                        end
                        NumPoints_C{loc} = NumPoints{loc}(currentFlies,:);
                    end
                    [meanDist,meanDist2,sampDist,sampDist2,Z1,aStruct,Dn,otherInfo]...
                        = SyntheticFlyDecoder(data,model,likely_high_state_sorted,...
                        NumPoints,durVec,{B_C},NumPoints,{c2cons{cl}(currentFlies)},options,0);
                    S{cl}.(l{i}){j,c}.meanDist = meanDist{1};
                    S{cl}.(l{i}){j,c}.meanDist2 = meanDist2{1};
                    S{cl}.(l{i}){j,c}.sampDist = sampDist{1};
                    S{cl}.(l{i}){j,c}.sampDist2 = sampDist2{1};
                    S{cl}.(l{i}){j,c}.Z1 = Z1{1};
                    S{cl}.(l{i}){j,c}.aStruct = aStruct{1};
                    S{cl}.(l{i}){j,c}.Dn = Dn{1};
                    S{cl}.(l{i}){j,c}.otherInfo = otherInfo;
                    close all
                end
            end
        end
    end
end
save(options.synthKMFile,'S','-v7.3')

end