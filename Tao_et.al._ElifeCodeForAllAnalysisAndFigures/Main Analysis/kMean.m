function [kMeans] = kMean(distance,scenario,options,plotting,saveDat)
%kMean Takes in the B matrix with high level state distributions for
%individual flies used for generating the synthetic flies and computes the
%kMeans clustering that places the flies into manimum number of clusters.

% Inputs:
%    syntheticMat: matrix containing the high level states used for
%    synthetic fly generation
%    distance: 'Euclidean'
%    resultfile: file to print to
%    scenario: ordering of scenarios and locations
%    plot: 1 for yes and 0 for no
%    saveDat: 1 for yes and 0 for no

% Output:
%    kMeans: A structure with fields:
%       Seed: 1xn structures containing kMeans for all possible starting
%       seeds
%       maxCluster: structure containing kMeans based on seed that
%       maximizes the number of total clusters per scenario
%       scenarios: same as input
%       empirical: B matrix used in clustering

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) Probability as a function of clusters
%       2.) Flies in each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% 2017, Bhandawat Lab

if (isempty(options.synthFile))
    disp('Need Synthetic file')
end
if (~exist('distance', 'var'))
    distance = 'euclidean';
end
resultfile = 'PDF_Files/test';
if (~exist('scenario', 'var'))
    scenario = {'b_i','d_i','diff_i';'b_o','d_o','diff_o'};
end
if (~exist('plotting', 'var'))
    plotting = 0;
end
load(options.synthFile,'B_new','flyClustNdx')

kMeansAll = cell(size(flyClustNdx));

for i = 1:length(flyClustNdx)
    kMeans = [];
    maxNdx = cell(2,3);
    maxCluster = zeros(2,3);
    
    currFlies = flyClustNdx{i};
    for ndx = 1:length(currFlies)
        kk = 1;
        for state = 1:3
            for loc = 1:2
                % remove empty HLS
                tmp = B_new{i}{loc,state}';
                tmp(:,sum(tmp,1)==0)=[];
                features = num2cell(tmp,2);
                
                kMeansTemp = ml_kmeans_aicbic(features, resultfile, distance, ndx, 1);

                diff_bic=-diff(kMeansTemp.bic);
                df = length(features)-1;
                t = sqrt(log(df+1)+diff_bic);
                t_val = real(t);
                P_val = 2*tcdf(-abs(t_val),df);
                perCond = (kMeansTemp.bic-min(kMeansTemp.bic))./min(kMeansTemp.bic)<0.05;
                optCluster = find(P_val<0.05 & perCond(1:end-1),1,'last')+1;
                
                if isempty(optCluster)
                    optCluster = 1;
                end

                kMeansTemp.optCluster = optCluster;
                kMeansTemp.P_val = [NaN P_val];
                kMeansTemp.t_val = [NaN t_val];
                kMeans.Seed{ndx}.data{loc,state} = kMeansTemp;

                if optCluster>maxCluster(loc,state)
                    maxNdx{loc,state} = ndx;
                    maxCluster(loc,state) = optCluster;
                elseif optCluster == maxCluster(loc,state)
                    maxNdx{loc,state} = [maxNdx{loc,state} ndx];
                end

                if plotting == 1
                    figure(1)
                    subplot(3,2,kk)
                    plot(kMeansTemp.kmin:kMeansTemp.kmax,[1 P_val]);xlim([kMeansTemp.kmin,kMeansTemp.kmax])
                    text(10,0.7,['Optimal ' num2str(optCluster) ' clusters'])
                    text(10,0.6,['Pk = ' num2str(P_val(optCluster-1))])
                    xlabel('Clusters');ylabel('P');title(scenario{loc,state})

                    figure(2)
                    subplot(3,2,kk)
                    for j = 1:optCluster
                        currentCluster = kMeansTemp.cluster{optCluster}{j};
                        plot(currentCluster,j*ones(1,length(currentCluster)),'ok');hold on
                    end
                    axis([0 35 0 optCluster+1])
                    xlabel('States');ylabel('Cluster');title(scenario{loc,state})
                    kk = kk+1;
                end
            end
        end
        disp(ndx)
    end

    for state = 1:3
        for loc = 1:2
            kMeans.maxCluster.data{loc,state} = kMeans.Seed{maxNdx{loc,state}(end)}.data{loc,state};
            kMeans.maxCluster.data{loc,state}.seed = maxNdx(loc,state);
        end
    end
    kMeans.scenario = {'b_i','d_i','diff_i';'b_o','d_o','diff_o'};
    kMeans.empirical = B_new;
    
    kMeansAll{i} = kMeans;
end
kMeans = kMeansAll;

if saveDat == 1
    save(options.KMEmp,'kMeans','-v7.3')
end

end