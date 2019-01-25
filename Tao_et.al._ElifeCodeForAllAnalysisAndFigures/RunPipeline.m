close all
clear
clc
currentFolder = pwd;
addpath(genpath([currentFolder '/Support Files']))
addpath(genpath([currentFolder '/Cluster Validation']))
addpath(genpath([currentFolder '/Block Diagonalization']))
addpath(genpath([currentFolder '/Main Analysis']))
addpath(genpath([currentFolder '/Other Analysis']))
addpath(genpath([currentFolder '/Visualization Files']))

disp('Running pipeline...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating Data Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.datFile = 'Data/Aug6_WTACV0_perp_para_sd2_var.mat';
options.synthFile = 'mat_Files/Synthetic_100it.mat';
options.KMEmp = 'mat_Files/kMeansEmpirical.mat';
options.KMGLMFile = 'mat_Files/KMeansGLM.mat';
options.synthKMFile = 'mat_Files/Synthetic_100itKMCluster.mat';
options.empKMFile = 'mat_Files/kMeansEmpiricalGLM.mat';
options.numSynth = 100;

disp('Processing')
% Always run first
[data,model,~,likely_high_state_sorted,~,~,~,clusters_to_consider,~,~,~] = ...
    Processing(options);

disp('GLMALL')
% Do GLM on all flies
[NumPoints,durVec,indDurVec,type,HLSDistRawOn,HLSDistRawOff,PCHLSComp,bs_all,bad_Fly_List,flyClustNdx]...
    = GLM_decode(data,model,likely_high_state_sorted,clusters_to_consider);

disp('GLMVis')
fig_title = [];
% Visualize GLM results (used to generate B matrix)
[B] = GLM_Visualization(type,HLSDistRawOn,HLSDistRawOff,PCHLSComp,bs_all,bad_Fly_List,flyClustNdx,fig_title);
close all

disp('Synth Flies')
% Create synthetic tracks for all flies
% Generates Synthetic_100it.mat
[~,~,~,~,~,~,~,~] = SyntheticFlyDecoder(data,model,...
    likely_high_state_sorted,NumPoints,durVec,B,NumPoints,flyClustNdx,options,1);

disp('kMeans')
% Create initial kMeans clustering based on empirical HLS distribution used
% to generate synthetic flies
% Generates kMeansEmpirical.mat
distance = 'euclidean';
scenario = {'b_i','d_i','diff_i';'b_o','d_o','diff_o'};
[kMeans] = kMean(distance,scenario,options,[],1);

disp('kMeansGLM')
% Improve KMeans clustering by using glm
% Generates KMeansGLM.mat
tic
[kMNew] = GLMClusterHandle(likely_high_state_sorted,flyClustNdx,options);
toc

disp('kMeansDecode')
% Create synthetic tracks for clusters of flies based on KMeans
% Generates Synthetic_100itKMCluster.mat
SyntheticClusterDecode(data,model,likely_high_state_sorted,kMeans,B,indDurVec,durVec,NumPoints,flyClustNdx,options)

for i = 1:length(clusters_to_consider)
    kMeans{i}.maxCluster = kMNew{i};
end
save(options.empKMFile,'kMeans');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except options
% list of visualization files
folder_PDF = 'PDF_Files'; fig_title = 'PDF_Files/VisualizationAll2';

% do processing again. Not necessary, if you don't delete variables
[data,model,S,likely_high_state_sorted,HHMM_by_State_sorted,HHMM_lowlevel_sorted,...
    likely_high_state_by_fly,clusters_to_consider,ndx,ndxLL,dur_sorted] = ...
    Processing(options);
[~,~,~,type,HLSDistRawOn,HLSDistRawOff,PCHLSComp,bs_all,bad_Fly_List,flyClustNdx] = ...
    GLM_decode(data,model,likely_high_state_sorted,clusters_to_consider);

disp('Plotting Visualization files...')
% main visualization files
Visualization(data,S,model,likely_high_state_sorted,HHMM_by_State_sorted,...
    likely_high_state_by_fly,clusters_to_consider,ndx,ndxLL,folder_PDF,fig_title);
close all
GLM_Visualization(type,HLSDistRawOn,HLSDistRawOff,PCHLSComp,bs_all,bad_Fly_List,flyClustNdx,fig_title);
close all
SyntheticVisualization(fig_title,options,flyClustNdx)
close all
plotXYPositionTracks(data,model,likely_high_state_sorted,clusters_to_consider,fig_title)
close all
load('mat_Files/KMeansGLM.mat','PCrawALL','PCHLSALL','Ind','kMNew');
for i = 1:length(clusters_to_consider)
    GLMCluster_Vis(PCrawALL{i},PCHLSALL{i},Ind{i},kMNew{i},clusters_to_consider(i),fig_title)
end

% other analysis/visualizations
% plot subsampling
tBins = [5 30:30:5400];
subSamplingHLSDist(data,model,likely_high_state_sorted,flyClustNdx,options,tBins,2,fig_title);
% plot time based analyses
TemporalAnalysis(data,fig_title);
SpaceTimeAnalysis(data,model,likely_high_state_sorted,flyClustNdx,fig_title);
% plotCurvature
CurvatureDistributions(model,HHMM_by_State_sorted,fig_title)
% Compare to random points in simplex
simplexXMeans(model,options,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
