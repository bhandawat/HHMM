function SyntheticVisualization(fig_title,options,flyClustNdx)
%SyntheticVisualisation takes in synthetic fly data matrices to generate
%visuals that interpret the results.

% Inputs:
%    SynthFName: the name of the data file containing the synthetic tracks
%    based on median empirical fly
%    kMfName: the name of the data file containing the empirical kMeans 
%    clustering
%    kMGLMfName: the name of the data file containing the synthetic tracks
%    based on kMeans clustering followed by GLM analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) visualizeclusters
%           a.) Plot 2D PC representation of GLM clusters after PCA
%       2.) HLSdistributions
%           a.) Mean HLS distributions for synthetic and for empirical flies
%       3.) plotDistances
%           a.) Analysis plots based on distances
%       4.) plotDistances2
%           a.) Compare empirical and synthetic histograms
%           b.) Box plot for median of synthetic tracks and empirical
%           c.) Correlation plot for distance between flies
%       5.) DendrogramVarianceAnalysis
%           a.) Ignore in curent version, update later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% 2017, Liangyu Tao

if (isempty(options.synthFile))
    disp('need synthetic file');
    return
end
if (isempty(options.empKMFile))
    disp('need kmeans GLM empirical file')
    return
end

if (isempty(options.synthKMFile))
    disp('need kmeans synthetic file')
    return
end

load(options.empKMFile,'kMeans')
load(options.synthFile,'B_new','Dn','aStruct','meanDist')

% plotting based on the synthetic tracks based on median empirical fly
% grouped into clusters
[Dn2] = plotDistancesKM(options,fig_title);close all

% plotting based on the synthetic tracks based on median empirical fly
flyClustNdx = flyClustNdx(~cellfun(@isempty,flyClustNdx));
Dn = Dn(~cellfun(@isempty,Dn));
aStruct = aStruct(~cellfun(@isempty,aStruct));
meanDist = meanDist(~cellfun(@isempty,meanDist));
for cl = 1:length(Dn)
    DnC = Dn{cl};aStructC = aStruct{cl};meanDistC = meanDist{cl};
    if ~isempty(DnC)
        plotDistances2(DnC,fig_title);close all
        visualizeclusters(B_new{cl},kMeans{cl},fig_title);close all
        visualizeIndividualFlyDists(B_new{cl},kMeans{cl},flyClustNdx{cl},fig_title);close all
        HLSdistributions(aStructC,meanDistC,DnC,fig_title);close all
    end
end

end

function [] = plotDistancesKMGLM(kMGLMfName,fig_title)
load(kMGLMfName)
l = {'Inside','Outside'};
numTrials = 100;
for loc = 1:2
    for cluster = 1:length(S.(l{loc}))
        close all
        Dn = S.(l{loc}){cluster}.Dn;
        aStruct = S.(l{loc}){cluster}.aStruct;
        meanDist = S.(l{loc}){cluster}.meanDist;
        meanDist2 = S.(l{loc}){cluster}.meanDist2;
        sampDist = S.(l{loc}){cluster}.sampDist;
        sampDist2 = S.(l{loc}){cluster}.sampDist2;
        B_new = S.(l{loc}){cluster}.otherInfo.B_new;
        leng = reshape(S.(l{loc}){cluster}.otherInfo.medLeng',1,4);
        leng_to_consider = 1:4;
%         
%         figTitle = ['T4/Synth_KMCluster_' (l{loc}) '_C' num2str(cluster)];
        
        plotDistances(Dn,aStruct,meanDist,meanDist2,numTrials,leng,leng_to_consider,fig_title)
        plotDistances2(Dn,fig_title)
    end
end

end

%%Fig 6 and Supplementary 6
function [Dn] = plotDistancesKM(options,fig_title)
load(options.synthFile,'B_new')
load(options.empKMFile,'kMeans')
load(options.synthKMFile,'S')
clearvars -except B_new kMeans S fig_title
S = S(~cellfun(@isempty,S));
l = {'Inside','Outside'};
NumTrials = 100;

for cl = 1:length(S)
    leng = reshape(S{cl}.Inside{1, 1}.otherInfo{1}.medLeng,1,4);
    for state = 1:2
        for loc = 1:2
            aStruct.synthTrackAll{state,loc} = [];
            optCluster = kMeans{cl}.maxCluster.data{loc,state}.optCluster;
            for c = 1:optCluster
                aStruct.synthTrackAll{state,loc} = [aStruct.synthTrackAll{state,loc};...
                    S{cl}.(l{loc}){state,c}.aStruct.synthTrackAll{state,loc}];
            end
            for nMed = 1:4
                for nIt = 1:NumTrials
                    aStruct.SynthDist{state,loc}{nMed,nIt} = [];
                    for c = 1:optCluster
                        aStruct.SynthDist{state,loc}{nMed,nIt} = ...
                            [aStruct.SynthDist{state,loc}{nMed,nIt};...
                            S{cl}.(l{loc}){state,c}.aStruct.SynthDist{state,loc}{nMed,nIt}];
                    end
                end
            end
        end
    end
    numFlies = size(B_new{cl}{1},2);
    Dn = calcDistances(aStruct,B_new{cl},NumTrials,leng,numFlies);
    plotDistances2(Dn,fig_title)
    close all
end
end

%%Fig 6 and Supplementary 6
function [] = visualizeclusters(B_new,kMeans,fig_title)
kk = 1;
colors = 'rkgcmrkgcmrkgcmrkgcmrkgcmrkgcm';
nHLS = size(B_new{1},1);
GlobalOptCluster = max([kMeans.maxCluster.data{1,1}.optCluster,...
    kMeans.maxCluster.data{1,2}.optCluster,...
    kMeans.maxCluster.data{2,1}.optCluster,...
    kMeans.maxCluster.data{2,2}.optCluster]);
for state = 1:2
    for loc = 1:2
        [~,score,~,~,~,~] = pca(B_new{loc,state}');
        optCluster=kMeans.maxCluster.data{loc,state}.optCluster;
        clabels = zeros(size(score,1),1);
        cluster = cell(1,optCluster);a = cell(1,optCluster);b = cell(1,optCluster);
        c = cell(1,optCluster);scores = cell(1,optCluster);raw = cell(1,optCluster);
        dx = 0.01; dy = 0.01;
        for i = 1:optCluster
            cluster{i}=kMeans.maxCluster.data{loc,state}.cluster{1,optCluster}{1,i};
            clabels(cluster{i}) = i;
            a{i}=[cluster{i}]';b{i}=num2str(a{i});c{i}=cellstr(b{i});
            scores{i}=score(cluster{i},:);
            raw{i} = B_new{loc,state}(:,cluster{i});
            figure(99);
            subplot(2,2,kk)
            text(scores{i}(:,1)+dx, scores{i}(:,2)+dy,c{i},'FontSize',8);hold on;
            
            figure(100);
            subplot(4,GlobalOptCluster,((state-1)*2+loc-1)*GlobalOptCluster+i)
            bar(mean(raw{i}'),'LineWidth',0.75)
            xlim([0 nHLS+1]);ylim([0 0.7])
            text(1.5,0.5,['Cluster ' num2str(i) ' (n=' num2str(length(c{i})) ')'],'Color',colors(i));hold off
            if i == 1
                title(kMeans.scenario(loc,state))
            end
            set(gcf,'position',[849 49 824 918])
        end
        figure(99);
        gscatter(score(:,1),score(:,2),clabels,colors(1:optCluster),'o')
        set(gcf,'position',[782 195 834 722])
        hold off;
        title(kMeans.scenario(loc,state))
        xlabel('Principal Component 1');
        ylabel('Principal Component 2');
        if loc == 1
            axis([-0.4 0.8 -0.4 0.3])
        else
            axis([-0.3 0.6 -0.3 0.3])
        end
        kk = kk+1;
    end
end
hold off;

GlobalOptCluster = max([kMeans.maxCluster.data{1,3}.optCluster,...
    kMeans.maxCluster.data{1,3}.optCluster]);
for state = 3
    for loc = 1:2
        optCluster=kMeans.maxCluster.data{loc,state}.optCluster;
        clabels = zeros(size(score,1),1);
        cluster = cell(1,optCluster);a = cell(1,optCluster);b = cell(1,optCluster);
        c = cell(1,optCluster);scores = cell(1,optCluster);raw = cell(1,optCluster);
        dx = 0.01; dy = 0.01;
        for i = 1:optCluster
            cluster{i}=kMeans.maxCluster.data{loc,state}.cluster{1,optCluster}{1,i};
            clabels(cluster{i}) = i;
            a{i}=[cluster{i}]';b{i}=num2str(a{i});c{i}=cellstr(b{i});
            scores{i}=score(cluster{i},:);
            raw{i} = B_new{loc,state}(:,cluster{i});
            
            figure(101);
            subplot(GlobalOptCluster,2,(i-1)*2+loc)
            bar(mean(raw{i}'),'LineWidth',0.75)
            xlim([0 nHLS+1]);ylim([-0.4 0.4]);
            text(3,0.3,['Cluster ' num2str(i) ' (n=' num2str(size(raw{i},2)) ')'],'Color',colors(i));hold off
            if i == 1
                title(kMeans.scenario(loc,state))
            end
            set(gcf,'position',[849 49 824 918])
        end
    end
end
hold off;
if ~isempty(fig_title)
    figure(99)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    figure(100)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    figure(101)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

end

function [] = visualizeIndividualFlyDists(B_new,kMeans,c2cons,fig_title)
nFlyPerPage = 9;
colors = 'rkgcmrkgcmrkgcmrkgcmrkgcmrkgcm';
nHLS = size(B_new{1},1);
for state = 1:3
    for loc = 1:2
        allDist = B_new{loc,state}';
        optCluster=kMeans.maxCluster.data{loc,state}.optCluster;
        clust = kMeans.maxCluster.data{loc,state}.cluster{optCluster};
        flyLable = cell(size(allDist,1),1);c = cell(size(allDist,1),1);
        for i = 1:optCluster
            flyLable(clust{i}) = {['Cluster ' num2str(i)]};
            c(clust{i}) = {colors(i)};
        end
        
        for fly = 1:size(allDist,1)
            rowNdx = mod(fly,nFlyPerPage);
            rowNdx(rowNdx==0) = nFlyPerPage;
            figNdx = ceil(fly./nFlyPerPage);
            
            figure(figNdx)
            subplot(nFlyPerPage,6,(rowNdx-1)*6+((state-1)*2+loc))
            bar(allDist(fly,:));hold on
            xlim([0 nHLS+1]);
            if state == 3
                ylim([-0.5 0.5]);
                text(3,0.3,flyLable{fly},'Color',c{fly});hold off
            else
                ylim([0 1]);
                text(3,0.75,flyLable{fly},'Color',c{fly});hold off
            end
            set(gcf,'position',[849 49 824 918])
            if rowNdx == 1
                title(kMeans.scenario(loc,state))
            end
            if ((state-1)*2+loc)==1
                ylabel(['Fly ' num2str(c2cons(fly))],'FontSize',11)
            end
        end
    end
end

if ~isempty(fig_title)
    for i = 1:figNdx
        figure(i)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
end

end


%%Supplementary 5
function [] = HLSdistributions(aStruct,meanDist,Dn,fig_title)
nHLS = size(meanDist{1},1);
for state = 1:2
    for loc = 1:2
        synth = zeros(1,nHLS);
        for i=1:nHLS
            synth(i)=length(find(aStruct.synthTrackAll{(state-1)*2+loc,1}{state,loc}==i));
        end
        synth = synth./(sum(synth));
        %synth = mean(aStruct.SynthDist{state,loc}{1,(state-1)*2+loc});
        
        figure(5)
        subplot(2,2,(loc-1)*2+state)
        bar(synth,'k');
        xlim([0 nHLS+1]);ylim([0 0.8]);
        title(Dn.scenario(state,loc))
        figure(6)
        subplot(2,2,(loc-1)*2+state)
        bar(meanDist{state,loc},'k')
        xlim([0 nHLS+1]);ylim([0 0.8]);
        title(Dn.scenario(state,loc))
    end
end
figure(5);suptitle('Synthetic');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end
figure(6);suptitle('Empirical');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

end

function [] = plotDistances(Dn,aStruct,meanDist,meanDist2,numTrials,leng,leng_to_consider,fig_title)
%close all
figure
b_i = leng_to_consider(1);
b_o = leng_to_consider(2);
d_i = leng_to_consider(3);
d_o = leng_to_consider(4);

% examine during inside scenario only
temp = zeros(1,10);
for i = 1:numTrials
    temp = temp+aStruct.meanDistNew{2,1,d_i}(i,:);                         % 4 is trial number
end
temp = temp/i;
subplot(2,1,1);bar(meanDist{2,1});xlim([0 11]);ylim([0 1])
title('Emperical During Inside');xlabel('HLS');ylabel('Probability')
subplot(2,1,2);bar(temp);xlim([0 11]);ylim([0 1])
title(['Synthetic During Inside ' num2str(leng(end)) ' pts']);xlabel('HLS');ylabel('Probability')
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

figure
ndxTemp = find(meanDist2{2,1,2}>0);
frameNdx = 1:1:length(meanDist2{2,1,2});
subplot(2,1,1);plot(frameNdx(ndxTemp),meanDist2{2,1,2}(ndxTemp),'.');
title('HLS 2');xlabel('Frames');ylabel('Occurances')
ndxTemp = find(meanDist2{2,1,3}>0);
frameNdx = 1:1:length(meanDist2{2,1,3});
subplot(2,1,2);plot(frameNdx(ndxTemp),meanDist2{2,1,3}(ndxTemp),'.');
title('HLS 3');xlabel('Frames');ylabel('Occurances')
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% plot Jensen-Shannon Divergence distribution for synthetic
figure
for l = leng_to_consider
    meanNjs = reshape(cell2mat(Dn.jsSyn.d_i.all(l,:)),length(Dn.jsSyn.d_i.all{l,1}),numTrials);
    meanNjs = sum(meanNjs,2)/numTrials;
    subplot(2,2,l);bar(meanNjs);xlim([0 21]);ylim([0 1])
    xlabel('Djs');ylabel('Probability')
    set(gca,'xticklabel',cellfun(@num2str, num2cell(0:1/4:1), 'UniformOutput', false))
    title(['Avg ' num2str(leng(l)) ' Frames'])
end
suptitle('Synthetic JS Divergence During Inside')
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% plot Euclidean Divergence distribution for synthetic
figure
for l = leng_to_consider
    meanNij = reshape(cell2mat(Dn.ijSyn.d_i.all(l,:)),length(Dn.ijSyn.d_i.all{l,1}),numTrials);
    meanNij = sum(meanNij,2)/numTrials;
    subplot(2,2,l);bar(meanNij);xlim([0 21]);ylim([0 1])
    xlabel('Dij');ylabel('Probability')
    set(gca,'xticklabel',cellfun(@num2str, num2cell(0:3/4:3), 'UniformOutput', false))
    title(['Avg ' num2str(leng(l)) ' Frames'])
end
suptitle('Synthetic Euclidean Distance During Inside')
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% plot Jensen-Shannon Divergence distribution for emperical
figure
kk = 1;
for state = 1:2                                             % 1 = before, 2 = during
    for loc = 1:2                                           % 1 = inside, 2 = outside
        meanNjs = Dn.jsEmp.all{state,loc};
        subplot(2,2,kk);bar(meanNjs);xlim([0 21]);ylim([0 1])
        xlabel('Djs');ylabel('Probability')
        set(gca,'xticklabel',cellfun(@num2str, num2cell(0:1/4:1), 'UniformOutput', false))
        title(Dn.scenario{state,loc})
        kk = kk+1;
    end
end
suptitle('Emperical JS Divergence')
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% plot Euclidean Divergence distribution for emperical
figure
kk = 1;
for state = 1:2                                             % 1 = before, 2 = during
    for loc = 1:2                                           % 1 = inside, 2 = outside
        meanNij = Dn.ijEmp.all{state,loc};
        subplot(2,2,kk);bar(meanNij);xlim([0 21]);ylim([0 1])
        xlabel('Djs');ylabel('Probability')
        set(gca,'xticklabel',cellfun(@num2str, num2cell(0:3/4:3), 'UniformOutput', false))
        title(Dn.scenario{state,loc})
        kk = kk+1;
    end
end
suptitle('Emperical Euclidean Distance')
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end


end

function [] = plotDistances2(Dn,fig_title)
% plot Euclidean Distance histogram for emperical and Euclidean Distance
% curve for synthetic
figure
kk = 1;
for state = 1:2                                             % 1 = before, 2 = during
    for loc = 1:2                                           % 1 = inside, 2 = outside
        SampRun_synth = [];
        numIt = length(Dn.ijSyn.(Dn.scenario{state,loc}).tot{(state-1)*2+loc});
        for it = 1:numIt
            SampRun_synth = [SampRun_synth; squareform(Dn.ijSyn.(Dn.scenario{state,loc}).tot{(state-1)*2+loc}{it})];% 3 is scenario, 2 is trial
        end
        Emp_All = squareform(Dn.ijEmp.tot{state,loc});
        %pval = ranksum(reshape(SampRun_synth,1,numel(SampRun_synth)),Emp_All);
        pval = ranksum(SampRun_synth(1,:),Emp_All);
        for i = 1:numIt
            p{state,loc}(i) = ranksum(SampRun_synth(i,:),Emp_All);
        end
        dx = 0.02;
        bins = [0:dx:1];
        y_synAll = zeros(numIt,length(bins));
        for it = 1:numIt
            Nsyn = histcounts(SampRun_synth(it,:),bins)/length(SampRun_synth(it,:));
            pd = fitdist(SampRun_synth(it,:)','Kernel');
            y_synAll(it,:) = pdf(pd,bins);
        end
        y_syn = mean(y_synAll,1);
        Nemp = histcounts(Emp_All,bins)/length(Emp_All);
        subplot(2,2,kk)
        bar(bins(2:end)-dx,Nemp,'FaceColor',[0.1 0.1 0.1]);hold on
        plot(bins,y_syn/sum(y_syn),'r','LineWidth',1)
        kk = kk+1;
        xlim([0 1]);ylim([0 dx*12.5])
        title([Dn.scenario{state,loc} ' ' num2str(pval)])
    end
end
%print('-dpsc2',['PDF_Files/Synth_decode/FigPanels' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

%plotting median distance box plots for synthetic and point for emperical
figure
kk = 1;
for state = 1:2                                             % 1 = before, 2 = during
    for loc = 1:2                                           % 1 = inside, 2 = outside
        subplot(2,2,kk);boxplot(cell2mat(Dn.ijSyn.(Dn.scenario{state,loc}).med((state-1)*2+loc,:)));hold on
        plot(Dn.ijEmp.med{state,loc},'bo')
        ylim([0 0.5])
        xlim([0.9 1.1])
        ylabel('Median Euclidean Distance')
        title(Dn.scenario{state,loc})
        kk = kk+1;
        [~,pVal]=ttest2(cell2mat(Dn.ijSyn.(Dn.scenario{state,loc}).med((state-1)*2+loc,:)),Dn.ijEmp.med{state,loc});
        %pVal = signrank(cell2mat(Dn.ijSyn.(Dn.scenario{state,loc}).med((state-1)*2+loc,:)),Dn.ijEmp.med{state,loc});
        text(1,0.2,num2str(pVal))
    end
end
%print('-dpsc2',['PDF_Files/Synth_decode/FigPanels' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

figure
for loc = 1:2                                           % 1 = inside, 2 = outside
    durEmpD = squareform(Dn.ijEmp.tot{2,loc});%during
    befEmpD = squareform(Dn.ijEmp.tot{1,loc});%before
    p = polyfit(befEmpD,durEmpD,1);
    yfit = polyval(p,[0 befEmpD 1]);
    
    yresid = durEmpD - yfit(2:end-1);
    SSresid = sum(yresid.^2);
    SStotal = (length(durEmpD)-1) * var(durEmpD);
    rsq = 1 - SSresid/SStotal;
    rsq_adj = 1 - SSresid/SStotal * (length(durEmpD)-1)/(length(durEmpD)-length(p));
    subplot(2,1,loc)
    plot(befEmpD,durEmpD,'bo','MarkerSize',4);hold on
    plot([0 befEmpD 1],yfit)
    xlabel(['Distance bet flies before'])
    ylabel('Distance bet flies during')
    
    [R,P] = corrcoef(befEmpD,durEmpD);
    if loc == 1
        title(['Inside Rsq=' num2str(rsq_adj)])
    else
        title(['Outside Rsq=' num2str(rsq_adj)])
    end
    xlim([0 1])
    ylim([0 1])
end
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/FigPanels' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

end

function [] = DendrogramVarianceAnalysis(aStruct,numTrials,leng,leng_to_consider,fig_title)
% Plotting the maximum distance of all simulations
kk = 1;
leg = cellfun(@num2str, num2cell(leng), 'UniformOutput', false);
meanMax = cell(2,2);
s = {'before', 'during'};
l = {'inside', 'outside'};
figure
for loc = 1:2
    for state = 1:2
        meanMax{state,loc} = mean(aStruct.maxDistance{state,loc},2);
        subplot(2,2,kk);hold on
        for i = leng_to_consider
            plot(aStruct.maxDistance{state,loc}(i,:))
        end
        lgd = legend(leg,'Location','NorthEast');
        title(lgd,'Num Samples')
        xlabel('Trial');ylabel('Max Distance');xlim([0 numTrials]);ylim([0 6])
        title([s{state} ' ' l{loc}])
        kk = kk+1;
    end
end
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end


% Plotting the maximum distance of all simulations
kk = 1;
leg = cellfun(@num2str, num2cell(leng), 'UniformOutput', false);
s = {'before', 'during'};
l = {'inside', 'outside'};
figure
for loc = 1:2
    for state = 1:2
        [len_sorted,len_idx] = sort(leng);
        subplot(2,2,kk);plot(len_sorted,meanMax{state,loc}(len_idx));hold on
        for i = 1:length(leng)
            plot(leng,meanMax{state,loc},'r.');
        end
        xlabel('Number of Samples');ylabel('Avg Max Distance');xlim([0 max(leng)]);ylim([0 6])
        title([s{state} ' ' l{loc}])
        kk = kk+1;
    end
end
set(gcf,'position',[9 49 800 918])
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% plotting the average distances in a table
meanMax = [meanMax{1,1},meanMax{1,2},meanMax{2,1},meanMax{2,2}];
Before_In = meanMax(:,1);
Before_Out = meanMax(:,2);
During_In = meanMax(:,3);
During_Out = meanMax(:,4);
T = table(Before_In(len_idx),Before_Out(len_idx),During_In(len_idx),During_Out(len_idx),'RowNames',leg(len_idx)');
figure;
% Get the table in string form.
TString = evalc('disp(T)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
%print('-dpsc2',['PDF_Files/Synth_decode/Synth_and_Emp_Dist' '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end
end
