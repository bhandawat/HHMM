function [] = GLMCluster_Vis(PCrawALL,PCHLSALL,Ind,kMeans,c2cons,fig_title)
%GLMCluster_Vis Takes in outputs from GLMCluster and provides visualization
%to compare GLM results of clusters against all flies

% Inputs:
%    PCrawALL: Cell containing the percent correct predictions of raw
%    observables for all flies
%    PCHLSALL: Cell containing the percent correct predictions of HLS for
%    all flies
%    Ind: Cell containing all GLM information for individual fly decoding
%    kMeans: Cell containing kMeans clustering of flies.
%    fname: file name of .ps file to print to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) Group GLM results for each cluster HLS
%       2.) Group GLM results for each cluster HLS vs Raw
%       3.) Box plots showing distance to line clusters separately
%       4.) Box plots showing distance to line all clusters
%       5.) GLM results for each individual HLS vs cluster HLS
%       6.) GLM results for each all flies HLS vs cluster HLS
%       7.) Box plots showing distance to line comparing All,Clusters,Ind
%       9.) Individual GLM results for HLS vs Raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% 2017, Bhandawat Lab

cc=jet(10);
close all
tit1 = {'Inside','Outside'};
for i = 1:length(c2cons)
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        figure(locK);
        kk = 1;
        for clusters = 1:numCluster
            subplot(5,2,kk)
            y = PCrawALL{c2cons(i),1}.KC{locK,clusters};
            x = 0.5*ones(1,length(y'));
            x = x+0.5*rand(size(x));
            plot(x,y','o','MarkerSize',4);hold on;
            
            y = PCrawALL{c2cons(i),2}.KC{locK,clusters};
            x = 0.5*ones(1,length(y'));
            x=x+0.5*rand(size(x))+1;
            plot(x,y','o','MarkerSize',4);hold on; plot([0.25 2.25], [0.5 0.5]);hold off
            ylim([0 1])
            xlim([0.4 2.1])
            %title(['Group Fit Raw ' tit1{locK} ' Cluster ' num2str(clusters)])
            title(['Group: Raw ' tit1{locK} ' GMM C: ' num2str(c2cons(i)) ' Diff, C:' num2str(clusters)])
            kk = kk+2;
        end
    end
end

tit1 = {'Inside','Outside'};
for i = 1:length(c2cons)
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        figure(locK);
        kk = 2;
        for clusters = 1:numCluster
            subplot(5,2,kk)
            y = PCHLSALL{c2cons(i),1}.KC{locK,clusters};
            x = 0.5*ones(1,length(y'));
            x = x+0.5*rand(size(x));
            plot(x,y','o','MarkerSize',4);hold on;
            
            y = PCHLSALL{c2cons(i),2}.KC{locK,clusters};
            x = 0.5*ones(1,length(y'));
            x=x+0.5*rand(size(x))+1;
            plot(x,y','o','MarkerSize',4);hold on; plot([0.25 2.25], [0.5 0.5]);hold off
            ylim([0 1])
            xlim([0.4 2.1])
            %title(['Group Fit HLS ' tit1{locK} ' Cluster ' num2str(clusters)])
            title(['Group: HLS ' tit1{locK} ' GMM C: ' num2str(c2cons(i)) ' Diff, C:' num2str(clusters)])
            kk = kk+2;
        end
        set(gcf,'position',[849 49 824 918])
        %print('-dpsc2',[fname '.ps'],'-loose','-append');
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end

% Plot Raw vs HLS predictions for group (7B)
for i = 1:length(c2cons)
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        figure(locK+2);
        kk = 1;
        for clusters = 1:numCluster
            for loc = 1:2
                subplot(5,2,kk)
                yRaw = PCrawALL{c2cons(i),loc}.KC{locK,clusters};
                yHLS = PCHLSALL{c2cons(i),loc}.KC{locK,clusters};
                scatter(yRaw',yHLS',10*ones(size(yHLS')),repmat(cc(1,:),length(yHLS'),1)); hold on
                ylim([0 1])
                xlim([0 1])
                title(['Group: ' tit1{locK} ' GMM C: ' num2str(c2cons(i)) ' Diff, C:' num2str(clusters) ' Loc: ' tit1{loc}])
                kk = kk+1;
                xlabel('Raw');ylabel('HLS')
                refline(1)
            end
        end
        set(gcf,'position',[849 49 824 918])
        %print('-dpsc2',[fname '.ps'],'-loose','-append');
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end


% Calculate distance to line predictions
p = cell(1,length(c2cons));pT = cell(1,length(c2cons));d2 = cell(1,length(c2cons));
for i = 1:length(c2cons)
    p{i} = cell(2,3);
    pT{i} = cell(2,1);
    d2{i} = cell(2,2);
    
    figure(i+4);
    kk = 1;
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        flyTot = 1;
        for clusters = 1:numCluster
            d = [];
            for loc = 1:2
                yRaw = PCrawALL{c2cons(i),loc}.KC{locK,clusters};
                yHLS = PCHLSALL{c2cons(i),loc}.KC{locK,clusters};
                for fly = 1:length(yRaw)
                    pt = [yRaw(fly),yHLS(fly),0];               % Raw,HLS
                    d(loc,fly) = dist_to_line(pt, [0 0 0], [1 1 0]);                % + = HLS is better
                    %d2{locK}(loc,flyTot) = dist_to_line(pt, [0 0 0], [1 1 0]);                % + = HLS is better
                    flyTot = flyTot+1;
                    d2{i}{locK,loc} = [d2{i}{locK,loc} d(loc,fly)];
                end
                [p{i}{locK,clusters}(loc),~,~] = signrank(d(loc,:)');
            end
            subplot(5,2,kk)
            boxplot(d',{'Inside','Outside'})
            text(0.75,0.8,num2str(p{i}{locK,clusters}(1)))
            text(1.75,0.8,num2str(p{i}{locK,clusters}(2)))
            ylim([-1,1])
            ylabel('distance')
            title(['Group: ' tit1{locK} ' GMM C: ' num2str(c2cons(i)) ' Diff, C:' num2str(clusters)])
            kk = kk+1;
        end
    end
    set(gcf,'position',[849 49 824 918])
    %print('-dpsc2',[fname '.ps'],'-loose','-append');
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
end

figure(i+5)
kk = 1;
for i = 1:length(c2cons)
    for locK = 1:2
        for loc = 1:2
            [pT{i}{locK}(loc),~,~] = signrank(d2{i}{locK,loc}');
        end
        subplot(2,1,kk)
        boxplot([d2{i}{locK,1};d2{i}{locK,1}]',{'Inside','Outside'})
        text(0.75,0.8,num2str(pT{i}{locK}(1)))
        text(1.75,0.8,num2str(pT{i}{locK}(2)))
        ylim([-1,1])
        ylabel('distance')
        kk = kk+1;
        title(['Group: ' tit1{locK} ' GMM C: ' num2str(c2cons(i)) ' Diff, C: All'])
    end
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
end

cc=varycolor(10);
PCraw = [];
PCHLS = [];
for i = 1:length(c2cons)
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        fly = 1;
        for clusters = 1:numCluster
            for f = 1:length(Ind{c2cons(i)}.KC{locK,clusters}.fly)
                for loc = 1:2
                    PCraw{i}{locK}{loc}(fly) = Ind{c2cons(i)}.KC{locK,clusters}.fly{f}.PCraw{loc};
                    PCHLS{i}{locK}{loc}(fly) = Ind{c2cons(i)}.KC{locK,clusters}.fly{f}.PC{loc};
                end
                fly = fly+1;
            end
        end
    end
end

% Ind vs clusters
for i = 1:length(c2cons)
    figure;
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        kk = 1;
        yINdx = 1;
        for clusters = 1:numCluster
            if locK == 1
                yC = PCHLSALL{c2cons(i),1}.KC{locK,clusters};
                %yINdx = kMeans.data{locK,3}.cluster{numCluster}{clusters};      %fly ndx in cluster 1
                yI = PCHLS{i}{locK}{1}(yINdx:yINdx+length(yC)-1);
                yINdx = yINdx+length(yC);
                subplot(2,1,1)
                scatter(yI',yC',10*ones(size(yI')),repmat(cc(clusters+1,:),length(yI'),1)); hold on
                refline(0,0.5);plot([0.5 0.5],[0 1])
                axis([0 1 0 1])
                refline(1)
                xlabel('PC Ind');ylabel('PC KMeans')
                title(['Inside kMeans Inside'])
            end
            if locK == 2
                yC = PCHLSALL{c2cons(i),2}.KC{locK,clusters};
                %yINdx = kMeans.data{locK,3}.cluster{numCluster}{clusters};      %fly ndx in cluster 1
                yI = PCHLS{i}{locK}{2}(yINdx:yINdx+length(yC)-1);
                yINdx = yINdx+length(yC);
                subplot(2,1,2)
                scatter(yI',yC',10*ones(size(yI')),repmat(cc(clusters+1,:),length(yI'),1)); hold on
                refline(0,0.5);plot([0.5 0.5],[0 1])
                axis([0 1 0 1])
                refline(1)
                hline.Color = 'r';
                xlabel('PC Ind');ylabel('PC KMeans')
                title(['Outside kMeans Outside'])
            end
        end
    end
    set(gcf,'position',[849 49 824 918])
    %print('-dpsc2',[fname '.ps'],'-loose','-append');
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
end

% All vs clusters
for i = 1:length(c2cons)
    figure;
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        kk = 1;
        yINdx = 1;
        for clusters = 1:numCluster
            if locK == 1
                yC = PCHLSALL{c2cons(i),1}.KC{locK,clusters};
                yA = PCHLSALL{c2cons(i),1}.AC{locK,clusters};
                yINdx = yINdx+length(yC);
                subplot(2,1,1)
                scatter(yA',yC',10*ones(size(yA')),repmat(cc(clusters,:),length(yA'),1)); hold on
                xlabel('PC All');ylabel('PC KMeans')
                title(['Inside kMeans Inside'])
            end
            if locK == 2
                yC = PCHLSALL{c2cons(i),2}.KC{locK,clusters};
                yA = PCHLSALL{c2cons(i),2}.AC{locK,clusters};
                yINdx = yINdx+length(yC);
                subplot(2,1,2)
                scatter(yA',yC',10*ones(size(yA')),repmat(cc(clusters,:),length(yA'),1)); hold on
                hline.Color = 'r';
                xlabel('PC All');ylabel('PC KMeans')
                title(['Outside kMeans Outside'])
            end
        end
        refline(0,0.5);plot([0.5 0.5],[0 1])
        axis([0.4 1 0.4 1])
        refline(1)
        if numCluster==4
            legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Location','northwest')
        else
            legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Location','northwest')
        end
    end
    set(gcf,'position',[849 49 824 918])
    %print('-dpsc2',[fname '.ps'],'-loose','-append');
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
end


% Calculate distance to line predictions
%pt = zeros(1,34);pt2 = zeros(1,34);pt3 = zeros(1,34);
dKA = cell(length(c2cons),1);
dKI = cell(length(c2cons),1);
dAI = cell(length(c2cons),1);
for i = 1:length(c2cons)
    kk = 1;
    for locK = 1:2
        numCluster = kMeans.data{locK,3}.optCluster;
        flyTot = 1;
        for clusters = 1:numCluster
            d = [];
            yINdx = 1;
            for loc = 1:2
                yHLSK = PCrawALL{c2cons(i),loc}.KC{locK,clusters};
                yHLSA = PCHLSALL{c2cons(i),loc}.AC{locK,clusters};
                yHLSI = PCHLS{i}{locK}{loc}(yINdx:yINdx+length(yHLSK)-1);
                yINdx = yINdx;
                for fly = 1:length(yHLSK)
                    pt = [yHLSK(fly),yHLSA(fly),0];               % Raw,HLS
                    dKA{i}{locK,loc}{clusters}(fly)=dist_to_line(pt, [0 0 0], [1 1 0]);                 % + = first is better
                    pt2 = [yHLSK(fly),yHLSI(fly),0];
                    dKI{i}{locK,loc}{clusters}(fly)=dist_to_line(pt2, [0 0 0], [1 1 0]);                % + = first is better
                    pt3 = [yHLSA(fly),yHLSI(fly),0];
                    dAI{i}{locK,loc}{clusters}(fly)=dist_to_line(pt3, [0 0 0], [1 1 0]);                % + = first is better
                    flyTot = flyTot+1;
                end
            end
            yINdx = yINdx+length(yHLSK);
        end
    end
end

p_AF = cell(1,2);
location = {'Inside','Outside'};
comp = {'All-KMeans','Ind-KMeans','Ind-All'};
for j = 1:length(c2cons)
    figure
    for loc = 1:2
        d_AF(1,:) = cell2mat(dKA{j}{loc,loc});
        d_AF(2,:) = cell2mat(dKI{j}{loc,loc});
        d_AF(3,:) = cell2mat(dAI{j}{loc,loc});
        for i = 1:3
            [p_AF{loc}(i),~,~] = signrank(d_AF(i,:));
            subplot(3,2,(i-1)*2+loc)
            boxplot(d_AF(i,:),location{loc})
            text(0.99,0.8,num2str(p_AF{loc}(i)))
            ylim([-1,1]);xlim([0.9 1.1])
            ylabel('distance')
            title(['Group: ' comp{i} ' GMM C: ' num2str(c2cons(j)) ' ' location{loc}])
        end
    end
    set(gcf,'position',[849 49 824 918])
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
    clearvars d_AF
end

sub = 1;
for n = 1:length(c2cons)
    figure
    for locK = 1:2
        if locK==5
            set(gcf,'position',[9 49 824 918])
            %print('-dpsc2',[['PDF_Files/GLM_SCvsHLS'] '.ps'],'-loose','-append');
            figure
            sub = 1;
        end
        
        yRaw = PCraw{n}{locK}{1}(:);
        yRaw(yRaw==0) = [];
        y = PCHLS{n}{locK}{1}(:);
        y(y==0) = [];
        %         subplot(4,2,sub*2-1);
        subplot(2,1,1)
        scatter(yRaw',y',10*ones(size(y')),repmat(cc(n,:),length(y'),1)); hold on
        xlabel('PC raw');ylabel('PC HLS')
        title(['In avg over: ' num2str(30)])
        axis([0.5 1 0.5 1])
        refline(1)
        yRaw = PCraw{n}{locK}{2}(:);
        yRaw(yRaw==0) = [];
        y = PCHLS{n}{locK}{2}(:);
        y(y==0) = [];
        %         subplot(4,2,sub*2);
        subplot(2,1,2)
        scatter(yRaw',y',10*ones(size(y')),repmat(cc(n,:),length(y'),1)); hold on
        xlabel('PC raw');ylabel('PC HLS2')
        title(['Out avg over: ' num2str(30)])
        axis([0.5 1 0.5 1])
        refline(1)
        sub = sub+1;
    end
    suptitle(['Individual Fit' ' GMM C: ' num2str(c2cons(n))])
end
set(gcf,'position',[849 49 824 918])
%print('-dpsc2',[fname '.ps'],'-loose','-append');
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end
end


