function [BAll] = GLM_Visualization(type,HLSDistRawOn,HLSDistRawOff,...
    PCHLSComp,bs_all,bad_Fly_List,c2cons,fig_title)
%GLM_Visualization Takes in outputs from GLM_decode and conduct dendrogram
%analysis and provide visualization to compare the GLM results to the
%empirical HLS distributions

% Inputs:
%    type: matrix containing GLM results
%    HLSDistRawOn: empirical HLS distributions for during odor
%    HLSDistRawOff: empirical HLS distributions for before odor
%    PCHLSComp: HLS distribution matrix for HLS dist after PCA
%    bs_all: array showing number of grouped points used for GLM
%    bad_Fly_List: array showing number of points considered for GLM
%    fig_title: filename that figures are written to

% Outputs:
%    B: 2x3 cell where rows are location (in,out), columns are for state
%    (During,Before,Difference). Each cell contains n HLS distributions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) Individual GLM fits for inside (blue) and outside (red)
%       2.) Individual GLM fits on HLS vs raw observables
%       3.) Distance from the line of unity and Wilcoxon signed rank test
%       4-6.) Same as above, but for grouped flies
%       7.) Average probability distributions for HLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that the output of this will go into other analysis
%
% 2017, Bhandawat Lab

cc=jet(10);
close all
for l = 1:2
    for n=1:length(type)
        if ~isempty(type{n})
            flyNdx = c2cons{n};
            for f=1:length(flyNdx)
                i = flyNdx(f);
                Chance{n}{l}(i)=0.5;
                if sum(i == bad_Fly_List)<1
                    for m = 1:length(bs_all)
                        PCraw{n}{m}{l}(i)=type{n}.fly{i}.PCraw{m}{l};
                        PCrawALL{n}{m}{l}(i)=type{n}.fly{i}.PCrawALL{m}{l};

                        PCHLS{n}{m}{l}(i)=type{n}.fly{i}.PCHLS{m}{l};
                        PCHLSALL{n}{m}{l}(i)=type{n}.fly{i}.PCHLSALL{m}{l};
                    end
                end
            end
        end
    end
end

% Plot HLS predictions for individual (Fig 7 A)
figure
for n = 1:length(type)
    if ~isempty(type{n})
        figure(n);
        for m = 1:length(bs_all)
            subplot(4,2,m)
            y = PCHLS{n}{m}{1}(:);
            y(y == 0) = [];
            x = 0.5*ones(1,length(y'));
            x = x+0.5*rand(size(x));
            plot(x,y','o','MarkerSize',4);hold on;
            y = PCHLS{n}{m}{2}(:);
            y(y == 0) = [];
            x = 0.5*ones(1,length(y'));
            x=x+0.5*rand(size(x))+1;
            plot(x,y','o','MarkerSize',4);hold on; plot([0.25 2.25], [0.5 0.5]);hold off
            ylim([0 1])
            xlim([0.4 2.1])
            title(['average over: ' num2str(bs_all(m)) ' points'])
            %print('-dpdf', ['PDF_Decode/Fig7A.pdf']);
        end
        set(gcf,'position',[9 49 824 918])
        suptitle(['Individual GLM Fit Cluster: ' num2str(n)])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end

% Plot Raw vs HLS predictions for individuals(7B)
sub = 1;
for n = 1:length(type)
    if ~isempty(type{n})
        figure(n+length(type));set(gcf,'position',[9 49 824 918])
        suptitle(['Individual GLM Fit Cluster: ' num2str(n)])
        for m = 1:length(bs_all)
            if m==5
                set(gcf,'position',[9 49 824 918])
                figure
                sub = 1;
            end
            yRaw = PCraw{n}{m}{1}(:);
            yRaw(yRaw==0) = [];
            y = PCHLS{n}{m}{1}(:);
            y(y==0) = [];
            subplot(2,1,1)
            scatter(yRaw',y',10*ones(size(y')),repmat(cc(n,:),length(y'),1)); hold on
            xlabel('PC raw');ylabel('PC HLS')
            title(['In avg over: ' num2str(bs_all(m))])
            axis([0.5 1 0.5 1])
            refline(1)
            yRaw = PCraw{n}{m}{2}(:);
            yRaw(yRaw==0) = [];
            y = PCHLS{n}{m}{2}(:);
            y(y==0) = [];
            subplot(2,1,2)
            scatter(yRaw',y',10*ones(size(y')),repmat(cc(n,:),length(y'),1)); hold on
            xlabel('PC raw');ylabel('PC HLS2')
            title(['Out avg over: ' num2str(bs_all(m))])
            axis([0.5 1 0.5 1])
            refline(1)
            sub = sub+1;
        end
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end


% Calculate distance to line predictions
p = cell(1,2);
for n = 1:length(type)
    if ~isempty(type{n})
        figure(2*length(type)+n)
        for m = 1:length(bs_all)
            d = zeros(34,2);
            for l = 1:2
                for i = 1:length(PCraw{n}{m}{l})
                    pt = [PCraw{n}{m}{l}(i),PCHLS{n}{m}{l}(i),0];               % Raw,HLS
                    d(i,l) = dist_to_line(pt, [0 0 0], [1 1 0]);                % + = HLS is better
                end
                %d(:,l) = d(:,l).*abs(d(:,l));
                [p{m}(l),~,~] = signrank(d(:,l));
                %[~,p{m}(l),~] = ttest(d(:,l));
            end
            %subplot(4,2,m)
            boxplot(d,{'Inside','Outside'})
            text(0.75,0.35,num2str(p{m}(1)))
            text(1.75,0.35,num2str(p{m}(2)))
            ylim([-0.1,0.4])
            ylabel('distance')
            title([num2str(bs_all(m)) ' pt avg'])
        end
        set(gcf,'position',[849 49 824 918])
        suptitle(['Individual GLM Fit Cluster: ' num2str(n)])
        %print('-dpdf', ['PDF_Files/Fig6/GLM_Wilcoxon_signed_rank_test2']);
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end


% Plot HLS predictions for All (Fig 7 A)
for n = 1:length(type)
    if ~isempty(type{n})
        figure(2*length(type)+n)
        for m = 1:length(bs_all)
            subplot(4,2,m)
            y = PCHLSALL{n}{m}{1}(:);
            y(y == 0) = [];
            x = 0.5*ones(1,length(y'));
            x = x+0.5*rand(size(x));
            plot(x,y','o','MarkerSize',4);hold on;
            y = PCHLSALL{n}{m}{2}(:);
            y(y == 0) = [];
            x = 0.5*ones(1,length(y'));
            x=x+0.5*rand(size(x))+1;
            plot(x,y','o','MarkerSize',4);hold on; plot([0.25 2.25], [0.5 0.5]);hold off
            ylim([0 1])
            xlim([0.4 2.1])
            title(['average over: ' num2str(bs_all(m)) ' points'])
            %print('-dpdf', ['PDF_Decode/Fig7A.pdf']);
        end
        set(gcf,'position',[9 49 824 918])
        suptitle(['Population GLM Fit Cluster: ' num2str(n)])
        %print('-dpdf', ['PDF_Files/GLM_HLSOneCluster']);
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end

% Plot Raw vs HLS predictions for All(7B)
figure
sub = 1;
for n = 1:length(type)
    if ~isempty(type{n})
        figure(3*length(type)+n)
        for m = 1:length(bs_all)
            if m==5
                set(gcf,'position',[9 49 824 918])
                %print('-dpsc2',[['PDF_Files/GLM_SCvsHLS'] '.ps'],'-loose','-append');
                figure
                sub = 1;
            end
            yRaw = PCrawALL{n}{m}{1}(:);
            yRaw(yRaw==0) = [];
            y = PCHLSALL{n}{m}{1}(:);
            y(y==0) = [];
            %         subplot(4,2,sub*2-1);
            subplot(2,1,1)
            scatter(yRaw',y',10*ones(size(y')),repmat(cc(n,:),length(y'),1)); hold on
            xlabel('PC raw');ylabel('PC HLS')
            title(['In avg over: ' num2str(bs_all(m))])
            axis([0 1 0 1])
            refline(1)
            yRaw = PCrawALL{n}{m}{2}(:);
            yRaw(yRaw==0) = [];
            y = PCHLSALL{n}{m}{2}(:);
            y(y==0) = [];
            %         subplot(4,2,sub*2);
            subplot(2,1,2)
            scatter(yRaw',y',10*ones(size(y')),repmat(cc(n,:),length(y'),1)); hold on
            xlabel('PC raw');ylabel('PC HLS2')
            title(['Out avg over: ' num2str(bs_all(m))])
            axis([0 1 0 1])
            refline(1)
            sub = sub+1;
        end
        set(gcf,'position',[9 49 824 918])
        suptitle(['Population GLM Fit Cluster: ' num2str(n)])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end


% Calculate distance to line predictions
figure
p = cell(1,2);
for n = 1:length(type)
    if ~isempty(type{n})
        for m = 1:length(bs_all)
            d = zeros(34,2);
            for l = 1:2
                for i = 1:length(PCraw{n}{m}{l})
                    pt = [PCrawALL{n}{m}{l}(i),PCHLSALL{n}{m}{l}(i),0];               % Raw,HLS
                    d(i,l) = dist_to_line(pt, [0 0 0], [1 1 0]);                % + = HLS is better
                end
                %d(:,l) = d(:,l).*abs(d(:,l));
                [p{m}(l),~,~] = signrank(d(:,l));
                %[~,p{m}(l),~] = ttest(d(:,l));
            end
            %subplot(4,2,m)
            boxplot(d,{'Inside','Outside'})
            text(0.75,0.35,num2str(p{m}(1)))
            text(1.75,0.35,num2str(p{m}(2)))
            ylim([-0.1,0.4])
            ylabel('distance')
            title([num2str(bs_all(m)) ' pt avg'])
        end
        set(gcf,'position',[849 49 824 918])
        suptitle(['Population GLM Fit Cluster: ' num2str(n)])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end

scenarios = {'during in','during out','before in','before out','diff in','diff out'};
BAll = cell(size(c2cons));
for i = 1:length(c2cons)
    if ~isempty(c2cons{i})
        flyNdx = c2cons{i};
        B = cell(2,3);
        for l = 1:2
            flyNdx(bad_Fly_List) = [];
            B1 = cat(2,HLSDistRawOn{flyNdx,l});
            B2 = cat(2,HLSDistRawOff{flyNdx,l});
            B3 = B1-B2;
            B(l,:) = {B1,B2,B3};
        end
        figure;
        k = 1;
        for state = 1:3
            for loc = 1:2
                subplot(4,2,k)
                bar(mean(B{loc,state},2))
                xlim([0 size(B{1},1)+1]);title(scenarios{k})
                if state ==3
                    ylim([-0.2 0.2])
                else
                    ylim([0 0.5])
                end
                k = k+1;
            end
        end
        
        for loc = 1:2
            subplot(4,2,k)
            bar(mean(B{loc,3},2)./mean(B{loc,2},2))
            xlim([0 size(B{1},1)+1]);ylim([-1 1]);
            title('Fractional Change')
            k = k+1;
        end
        set(gcf,'position',[849 49 824 918])
        suptitle(['Population GLM Fit Cluster: ' num2str(i)])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
        BAll{i} = B;
    end
end

end
