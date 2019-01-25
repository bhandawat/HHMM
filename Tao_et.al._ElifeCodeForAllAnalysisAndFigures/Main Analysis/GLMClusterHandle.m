function [kMNew] = GLMClusterHandle(likely_high_state_sorted,c2cons,options)
%GLMClusterHandle Improves clusterin from KMeans by using GLM. This is he
%handle function. It calls GLMClusterAnalysis, GLMCluster, and
%GLMCluster_Vis

% Inputs:
%    kMfName: the name of the data file containing the empirical kMeans
%    clustering
%    data_file: name of HHMM model/data file
%    likely_high_state_sorted: Sorted HLS tracks from processing
%    vis: 1 = generate figures from GLMCluster_Vis, 0 = no

%
% 2017, Liangyu Tao

options.datFile = 'Data/Aug6_WTACV0_perp_para_sd2_var.mat';
options.KMEmp = 'mat_Files/kMeansEmpirical.mat';

if (isempty(options.KMEmp))
    disp('Need KMEmp file')
else
    load(options.KMEmp,'kMeans')
end

if (isempty(options.datFile))
    disp('Need HHMM Data file')
else
    load(options.datFile,'data','model')
end

data_old = data;
data = post_hoc_preprocess_LT(data_old);

allCNdx = find(~cellfun(@isempty,c2cons));
kMNew = cell(1,length(kMeans));PCrawALL = cell(1,length(kMeans));
PCHLSALL = cell(1,length(kMeans));Ind = cell(1,length(kMeans));
for c = 1:length(kMeans)
    c2consTmp = cell(size(c2cons));
    c2consTmp(allCNdx(c)) = c2cons(allCNdx(c));
    kMeansC = kMeans{c}.maxCluster;
    [data,model,likely_high_state_sorted,kMNew{c}] = GLMClusterAnalysis(kMeansC,data,model,likely_high_state_sorted,c2consTmp);
    [~,~,PCrawALL{c},PCHLSALL{c},~,Ind{c}] = GLMCluster(data,model,likely_high_state_sorted,kMNew{c},c2consTmp);
end
save(options.KMGLMFile,'kMNew','PCrawALL','PCHLSALL','Ind','kMNew','-v7.3');

end


