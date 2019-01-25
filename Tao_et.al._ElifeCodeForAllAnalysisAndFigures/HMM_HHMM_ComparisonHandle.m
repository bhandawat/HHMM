close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HHMM and HMM comparison. Note that currently code only works if HMM only
% 1 GMM cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.datFile = 'Data/Aug6_WTACV0_perp_para_sd2_var.mat';
options.datFileHMM = 'Data/Sept10_WTACV0_HMM_50_var2.mat';
options.synthFile = 'mat_Files/Synthetic_100it.mat';
options.KMEmp = 'mat_Files/kMeansEmpirical.mat';
options.KMGLMFile = 'mat_Files/KMeansGLM.mat';
options.synthKMFile = 'mat_Files/Synthetic_100itKMCluster.mat';
options.empKMFile = 'mat_Files/kMeansEmpiricalGLM.mat';
options.numSynth = 100;
fig_title = 'PDF_Files/HHMM_HMM_Comparison';

% parameters for clustering based on information bottleneck
param.K = 10;
param.maxIt = 100;
param.tau = 10;
param.beta = 100;

% calculate the duration distribution for HHMM
[~,~,~,~,~,HHMM_lowlevel_sorted,~,~,~,~,dur_HHMM_sorted] = Processing(options);

%conduct processing on HMM data
[model,data,HMMsortNdx,likely_state_by_fly_sorted,dur_HMM_States] = ProcessingHMM(options);

% conduct clustering based on information bottleneck
dur_HMM_sorted = HMMBlockClustering(model,data,HMMsortNdx,likely_state_by_fly_sorted,param,[]);

% compare the duration distributions of states of interest
statesOfInterest = {[1:10],5,4};
durationPlots(dur_HHMM_sorted,dur_HMM_States,dur_HMM_sorted,statesOfInterest,HHMM_lowlevel_sorted,[])
%durationPlots(dur_HHMM_sorted,dur_HMM_States,dur_HMM_sorted,statesOfInterest,fig_title)

