function Visualization(data,S,model,likely_high_state_sorted,HHMM_by_State_sorted,likely_high_state_by_fly,clusters_to_consider,ndx,ndxLL,folder_PDF,fig_title)
%Visualization Takes in the processed information from Processing.m to
%generate figures relating to the visualization of data.

% Inputs:
%    data: data files
%    S: structure including comprehensive information about fly tracks
%    model: HHMM class containing model information
%    likely_high_state_by_fly: Unsorted HLS tracks (34x10800)
%    HHMM_by_State_sorted: Sorted HLS data in the form of separate tracks
%    clusters_to_consider: clusters in HHMM
%    ndx,ndxLL: sorting ndx based on speed/curvature
%    cc: colormap used for plotting
%    folter_PDF: folder to plot to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) Transpose Property Matrices
%       2.) Speed and curvature contour maps for HLS
%       3.) forward trajectories
%       4.) HLS state densities over circular and rectangular grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Outputs:
%    N/A

% Note that all inputs are generated from GLM_Decode. Before and during
% indexes are switched in this analysis as we designate before then during
%
% 2017, Siddhi Ozarkar, Liangyu Tao

% fig_title = [];folder_PDF = [];
close all
% supplementary figure 2
TransposeProperties(model,clusters_to_consider,ndx,ndxLL,fig_title)

close all
% figure 2a,3b, supplementary 3
[~,pLL,X,Y] = HLSObsModContour(data,model,likely_high_state_by_fly,clusters_to_consider,ndx,fig_title);
% close all
% CommonPointsAnalysis(commTracks,data,fig_title);
% load('tmp.mat')

close all
% figure 2b/c
states_of_interest = model.Qdim;len = 50;
forwardTrackSampleState(model,S,HHMM_by_State_sorted,states_of_interest,len,ndxLL,ndx,pLL,X,Y,fig_title)

close all
% Figure 3a
forwardTrack(model,S,HHMM_by_State_sorted,folder_PDF)

close all
% Figure 4 a/b
[before_norm,dur_norm,~] = HLSDensitySpace(data,model,likely_high_state_sorted,clusters_to_consider,fig_title);
close all
HLSDensitySpaceCircle(model,before_norm,dur_norm,clusters_to_consider,fig_title)


end
% Main Figures
function [emp_sorted,pLLAll,XAll,YAll] = HLSObsModContour(data,model,likely_high_state_by_fly,clusters_to_consider,newHLSndx,fig_title)
%%
conf = 0.85; % confidence bound for computing contour for empirical data from the model.
[~,cidx]=max(model.p);
obsdims = [];
for i = 1:length(model.obsTypes)
    obsdims=[obsdims model.obsTypes{i}.idx 17];
end

for i=1:model.NC
    idx=find(cidx==i);
    if(length(idx)>=10)
        m1=model.HHMMs{i};
        Qbar=zeros(m1.Qdim,1);
        for k=1:m1.Qdim
            Abar{k}=zeros(m1.dim,1);
        end
        for k=1:m1.Qdim
            emp{i,k}=[];
        end
        for j=idx
            dtemp=data{j}(obsdims,:);
            dtemp(end+1,:) = j;
            stemp = likely_high_state_by_fly(j,:);
            Qbar=Qbar+sum(m1.Qp{j},2);                                          % % chance for each HLS over all tracks 10x1 mat
            sltemp=cell(1,m1.Qdim);
            for k=1:m1.Qdim
                Abar{k}=Abar{k}+sum(m1.p{j}(m1.Aidx{k},:),2);
                emp{i,k}=[emp{i,k},dtemp(:,(stemp==k))];
                sltempAll = find(stemp==k);
                
                for kk = 1:length(sltempAll)
                    [~,temp{k}(:,kk)] = max(m1.p{j}((k-1)*m1.dim+1:k*m1.dim,sltempAll(kk)));
                    sltemp{k} = [sltemp{k} temp{k}(:,kk)];
                end
                
                for kkk = 1:m1.dim
                    empLL{i,k,kkk}=[];
                    empLL{i,k,kkk}=[empLL{i,k,kkk},dtemp(:,sltempAll(sltemp{k}==kkk))];
                end
            end
        end
        emp_sorted(i,:) = emp(i,newHLSndx{clusters_to_consider==i});
        Qbar=Qbar/sum(Qbar);                                                    % distribution of HLS for entire track
        for k=1:m1.Qdim
            Abar{k}=Abar{k}/sum(Abar{k});                                       % Distribution of lower level states
        end
        
        figure(1)
        kk=0;
        pLL = cell(m1.Qdim,m1.dim);
        pHL = cell(m1.Qdim,1);
        pHLS = cell(m1.Qdim,1);
        vHLS = cell(m1.Qdim,1);
        X = cell(1,model.Qdim);Y = cell(1,model.Qdim);
        for k=newHLSndx{clusters_to_consider==i}
            if(Qbar(k)<0.0001)
            else
                kk=kk+1;
                xmin=Inf;xmax=-Inf;
                ymin=Inf;ymax=-Inf;
                for n=1:m1.dim
                    mu{n}=m1.obsModels{n,k}.dists{1}.mu;
                    Sigma{n}=m1.obsModels{n,k}.dists{1}.ESigma;
                    invSigma{n}=inv(Sigma{n});
                    if(Abar{k}(n)<0.02)
                    else
                        xmin=min(xmin,mu{n}(1)-4*sqrt(Sigma{n}(1,1)));
                        xmax=max(xmax,mu{n}(1)+4*sqrt(Sigma{n}(1,1)));
                        ymin=min(ymin,mu{n}(2)-4*sqrt(Sigma{n}(2,2)));
                        ymax=max(ymax,mu{n}(2)+4*sqrt(Sigma{n}(2,2)));
                    end
                end
                
                %             X{k}=xmin:0.01:xmax;
                %             Y{k}=ymin:0.01:ymax;
                X{k}=-0.2:0.01:6;
                Y{k}=-2:0.01:2;
                
                x=X{k}'*ones(1,length(Y{k}));
                y=ones(length(X{k}),1)*Y{k};
                p=zeros(size(x));
                
                for n=1:m1.dim
                    lnp = -1/2*( (x-mu{n}(1)).^2*invSigma{n}(1,1) + (y-mu{n}(2)).^2*invSigma{n}(2,2) + 2*(x-mu{n}(1)).*(y-mu{n}(2))*invSigma{n}(1,2));
                    lnp = lnp-1/2*log(det(2*pi*Sigma{n}));
                    pLL{k,n} = exp(lnp);
                    p=p+Abar{k}(n)*exp(lnp);
                end
                
                idx2=randi(size(emp{i,k},2),1,min(1000,size(emp{i,k},2)));
                subplot(4,3,kk), scatter(emp{i,k}(1,idx2),emp{i,k}(2,idx2),1,'MarkerEdgeColor',[0.5 0.5 0.5]), hold on
                text(2,1.3,num2str(length(idx2)));
                p=p/sum(p(:));                                                  % normalize
                pHL{k} = p;
                m=sort(p(:),'descend');
                v=max(find(cumsum(m)<conf));
                v=m(v);
                pHLS{kk} = p;vHLS{kk} = v;
                subplot(4,3,kk), contour(X{k},Y{k},p',[v v],'LineWidth',2,'Color','k'), title(['HL State ',num2str(k),': ',num2str(round(Qbar(k)*100)),'%'])
                axis([-1,5,-1.5,1.5])
                hold off
                grid ON;
                drawnow
            end
        end
        suptitle(['Cluster=' num2str(i) ', 85% contour'])
        set(gcf,'position',[9 49 824 918])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
        pLLAll{i} = pLL;XAll{i} = X;YAll{i} = Y;
        
        figure(2)
        kk=0;
        for k=newHLSndx{clusters_to_consider==i}
            for n = 1:m1.dim
                p2 = pLL{k,n};
                p2=p2/sum(p2(:));                                                  % normalize
                m=sort(p2(:),'descend');
                %            v=max(find(cumsum(m)<0.6827));
                v=max(find(cumsum(m)<conf));
                v=m(v);
                
                kk=kk+1;
                if kk == 5*model.dim+1
                    set(gcf,'position',[9 49 824 918])
                    if ~isempty(fig_title)
                        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
                    end
                    figure
                    kk = 1;
                end
                
                subplot(5,model.dim,kk);scatter(empLL{i,k,n}(1,:),empLL{i,k,n}(2,:),1,'MarkerEdgeColor',[0.5 0.5 0.5]), hold on
                subplot(5,model.dim,kk); contour(X{k},Y{k},p2',[v v],'LineWidth',2,'Color','k');
                title(['LL State ',num2str(n),': ',num2str(round(Abar{k}(n)*100)),'%'])
                refline(0);
                grid ON;
                axis([-1,5,-1.5,1.5])
            end
        end
        set(gcf,'position',[9 49 824 918])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
        
        Xedges = [X{k}(1)-0.005:0.01:X{k}(end)+0.005];
        Yedges = [Y{k}(1)-0.005:0.01:Y{k}(end)+0.005];
        
        clearvars -except data model likely_high_state_by_fly clusters_to_consider...
            newHLSndx fig_title obsdims cidx i emp_sorted pLLAll XAll YAll conf
    end
    close all
end
end
function [] = forwardTrackSampleState(model,S,HHMM_by_State_sorted,states_of_interest,len,ndxLL,ndx,pLL,X,Y,fig_title)
%%
%Plot 30 random tracks of state 8 which are greater than 30 timesteps
%Plot the same 30 tracks by performing a transformation such that the
%initial velocity is oriented towards y axis
goodclusters = true(1,size(ndxLL,1));%~all(cellfun(@isempty,squeeze(HHMM_by_State_sorted(:,:,1,1))'));

cc=varycolor(model.dim);
if size(cc,1)==5
    cc(1,:)=[0 0 1];
    cc(2,:)=[0 0 0];
    cc(3,:)=[1 0.85 0];
    cc(4,:)=[1 0 0];
    cc(5,:)=[0 1 0];
end

fig_count=1;
figure(fig_count)
for k=1:size(HHMM_by_State_sorted,1)
    if goodclusters(k)
        for kk = 1:length(states_of_interest)
            allDuration = cellfun(@(C) size(C,2), S.LL.allTrack{k}{states_of_interest(kk)});
            longDuration = find(allDuration>len);                                     %chose random tracks greater than 50 timesteps
            currNdxLL = squeeze(ndxLL{k,states_of_interest(kk)});
            for i = 1:length(currNdxLL)
                cc_low_level(i,:) = cc(currNdxLL(i),:);
            end
            for kkk =  1:min(length(longDuration),20)
                longLL_temp=S.LL.allTrack{k}{states_of_interest(kk)}{longDuration(kkk)};
                brkpt_Long = [1 find(diff(longLL_temp) ~=0)];
                longTrack = S.data{k}{states_of_interest(kk)}{longDuration(kkk)};
                figure(fig_count)
                j=rem(kkk,9);
                if j==0
                    j=9;
                    fig_count=fig_count+1;
                end
                subplot(3,3,j);
                hold on
                if length(brkpt_Long)>1
                    for kkkk = 1:length(brkpt_Long)-1
                        plot(longTrack(1,brkpt_Long(kkkk):brkpt_Long(kkkk+1)),...
                            longTrack(2,brkpt_Long(kkkk):brkpt_Long(kkkk+1)),...
                            'Color',cc_low_level(longLL_temp(brkpt_Long(kkkk)+1),:))
                    end
                else
                    plot(longTrack(1,:),longTrack(2,:),'Color',cc_low_level(longLL_temp(1),:))
                end
                plot(longTrack(1,1),longTrack(2,1),'b*')
                axis([-1 1 -1 1])
                hold off
                if j == 9
                    set(gcf,'position',[849 100 824 824])
                    if ~isempty(fig_title)
                        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
                    end
                end
            end
        end
        set(gcf,'position',[849 100 824 824])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
        close all
    end
end


fig_count = fig_count+1;
for k= 1:size(HHMM_by_State_sorted,1)
    if goodclusters(k)
        for kk = 1:length(states_of_interest)
            allDuration = cellfun(@(C) size(C,2), S.LL.allTrack{k}{states_of_interest(kk)});
            longDuration = find(allDuration>len);
            newX = [];figure(fig_count)
            for kkk = 1:min(length(longDuration),20)
                longLL_temp=S.LL.allTrack{k}{states_of_interest(kk)}{longDuration(kkk)};
                brkpt_Long = [1 find(diff(longLL_temp) ~=0)];
                longTrack = S.data{k}{states_of_interest(kk)}{longDuration(kkk)};

                figure(fig_count)
                hold on
                longTrack(1,:) = longTrack(1,:)-longTrack(1,1);                                                     % make everything relative to the first point of trajectory
                longTrack(2,:) = longTrack(2,:)-longTrack(2,1);                                                     % (y)
                xy  = longTrack(1,:)+1i*longTrack(2,:);                                                         % change to polar coordinates for easier calculation
                vel_x=diff(longTrack(1,:)); vel_y=diff(longTrack(2,:));
                vel=vel_x+vel_y*1i;
                theta=angle(sum(vel(1:10)));                                                                     % assume that the first 8 steps define the curvature (i.e 1 and 8 are both in the direct forward direction)
                xy2=xy.*exp(-theta*1i+(pi/2)*1i);                                                               % shift by 90 degrees to point up instead of to the right
                newY = imag(xy2);
                newX=real(xy2);
                if length(brkpt_Long)>1
                    for kkkk = 1:length(brkpt_Long)-1
                        plot(newX(1,brkpt_Long(kkkk):brkpt_Long(kkkk+1)),...
                            newY(1,brkpt_Long(kkkk):brkpt_Long(kkkk+1)),...
                            'Color',cc_low_level(longLL_temp(brkpt_Long(kkkk)+1),:))
                    end
                else
                    if(length(newY(newY>0))>length(newY(newY<0))) && (max(abs(sqrt(vel_x.^2+vel_y.^2)))<0.7)
                        plot(real(xy2),imag(xy2),'Color',cc_low_level(longLL_temp(1),:));
                    end
                end
            end
            if ~isempty(newX)
                plot(newX(1,1),newY(1,1),'b*')
                ylim([-0.1 1]);
                xlim([-1 1]);
                set(gcf,'position',[9 220 1000 600])
            end
        end
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
        close all
    end
end

pLL(cellfun('isempty', pLL)) = [];
X(cellfun('isempty', X)) = [];
Y(cellfun('isempty', Y)) = [];
for i= 1:size(HHMM_by_State_sorted,1)
    kk = 1;figure;set(gcf,'position',[4 806 800 142]);
    if goodclusters(i)
        for ii = 1:length(states_of_interest)
            for k=ndx{i}(states_of_interest)
                for n = ndxLL{i,states_of_interest(ii)}

                    p2 = pLL{i}{k,n};
                    p2=p2/sum(p2(:));                                                  % normalize
                    m=sort(p2(:),'descend');
                    conf = 0.85;
                    v=max(find(cumsum(m)<conf));
                    v=m(v);
                    subplot(length(states_of_interest),5,kk); contour(X{i}{k},Y{i}{k},p2',[v v],'LineWidth',2,'Color',cc_low_level(n,:));
                    kk=kk+1;
                    refline(0);
                    grid ON;
                    axis([-1,5,-2,2])
                end
            end
            title(['State ' num2str(states_of_interest(ii))])
        end
        suptitle(['Cluster ' num2str(i)])
        if ~isempty(fig_title)
            print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        end
    end
end


end
function [] = forwardTrack(model,S,HHMM_by_State_sorted,folder_PDF)
%%
%this code takes the tracks of all states and plots them by performing a
%transformation such that their initial velocity is oriented toward y axis
%It then plots the median track for a particular state
if exist('fig_count','var')==0
    fig_count=1;
else
    fig_count=fig_count+1;
end
ColorSet=varycolor(100);

for k=1:size(HHMM_by_State_sorted,1)
    figure(fig_count)
    states_of_interest = 1:1:length(S.data{k});
    for kk=1:length(states_of_interest)
        h1=subplot(2,ceil(model.Qdim./2),kk);
        set(h1, 'ColorOrder', ColorSet);
        hold all;
        for kkk=1:length(S.data{k}{kk})
            x = S.data{k}{kk}{kkk}(1,:);
            y = S.data{k}{kk}{kkk}(2,:);
            x = x-x(11);                                                     % make everything relative to the first point of trajectory (x)
            y = y-y(11);                                                     % (y)
            xy  = x+1i*y;                                                   % change to polar coordinates for easier calculation
            vel_x=diff(x); vel_y=diff(y);
            vel=vel_x+vel_y*1i;
            theta=angle(sum(vel(1:10)));                                     % assume that the first 10 steps define the curvature (i.e 1 and 8 are both in the direct forward direction)
            xy2=xy.*exp(-theta*1i+(pi/2)*1i);                               % shift by 90 degrees to point up instead of to the right
            newY = imag(xy2);
            if(length(newY(newY>0))>length(newY(newY<0))) && (max(abs(sqrt(vel_x.^2+vel_y.^2)))<0.7)
                plot(real(xy2(11:end)),imag(xy2(11:end)));                                  % only plot tracks where the data doesn't jump around too much
            end
        end
        ylim([-0.1 1]);
        xlim([-1 1]);
        yticks([0 0.5 1])
        xticks([-1 0 1])
        title(['State ' num2str(states_of_interest(kk))])
    end
    set(gcf,'position',[9 220 1600 600])
    if ~isempty(folder_PDF)
        print('-painters','-depsc',[folder_PDF '/Fig3A_eps'])
    end
    fig_count=fig_count+1;
    
end
end
function [before_norm,during_norm,diff_norm] = HLSDensitySpace(data,model,likely_high_state_sorted,clusters_to_consider,fig_title)
%% plotting density of states as a function of space
if exist('fig_count','var')==0
    fig_count=1;
else
    fig_count=fig_count+1;
end

before_norm = cell(1,length(clusters_to_consider));
during_norm = cell(1,length(clusters_to_consider));
diff_norm = cell(1,length(clusters_to_consider));
for k=1:length(clusters_to_consider)
    cluster=clusters_to_consider(k);
    fly_no=find(model.p(cluster,:)>0.5);
    number_of_flies=length(fly_no);
    %bins=-1:0.02:1;
    bins=-1:1/30:1;
    before_density=zeros(length(bins),length(bins));
    before_state=zeros(length(bins),length(bins),model.Qdim);
    during_density=zeros(length(bins),length(bins));
    during_state=zeros(length(bins),length(bins),model.Qdim);
    for i=1:number_of_flies
        idx4=find(data{1,fly_no(i)}(11,:)==1 & data{1,fly_no(i)}(12,:)==1); % odor on and fly inside
        if isempty(idx4)==0
            first_entry = idx4(1);
        else
            first_entry = size(data{fly_no(i)},2);
        end
        % setting points outside of arena to at the arena boundary
        x=data{1,fly_no(i)}(1,:);
        x(x>1)=1;
        x(x<-1)=-1;
        y=data{1,fly_no(i)}(2,:);
        y(y>1)=1;
        y(y<-1)=-1;
        xbin=discretize(x,bins);
        ybin=discretize(y,bins);
        
        % calculating P(state|x,y) and P(all states|x,y)
        state=likely_high_state_sorted(fly_no(i),:);
        for ii=1:first_entry
            if state(ii)>0
                before_state(xbin(ii),ybin(ii),state(ii))=before_state(xbin(ii),ybin(ii),state(ii))+1;
                before_density(xbin(ii),ybin(ii))=before_density(xbin(ii),ybin(ii))+1;
            end
            
        end
        for iii=first_entry+1:length(xbin)
            if state(iii)>0
                during_state(xbin(iii),ybin(iii),state(iii))=1+during_state(xbin(iii),ybin(iii),state(iii));
                during_density(xbin(iii),ybin(iii))=during_density(xbin(iii),ybin(iii))+1;
            end
            
        end
    end
    figure(fig_count)
    before_state_normalized=NaN(length(bins),length(bins),model.Qdim);
    during_state_normalized=NaN(length(bins),length(bins),model.Qdim);
    for i=1:model.Qdim
        before_state_normalized(:,:,i)=rdivide(before_state(:,:,i),before_density); % Normalizing by doing P(state|x,y)/P(all states|x,y)
        before_state_normalized(isnan(before_state_normalized)) = 0 ;
        subplot(model.Qdim/2,4,(2*i-1))
        imagesc(before_state_normalized(:,:,i),[0 1]);
        caxis([-1 1]);
        
        during_state_normalized(:,:,i)=rdivide(during_state(:,:,i),during_density); % same, but after first entry
        during_state_normalized(isnan(during_state_normalized)) = 0 ;
        subplot(model.Qdim/2,4,2*i)
        imagesc(during_state_normalized(:,:,i),[0 1]);
        caxis([-1 1]);
    end
    set(gcf,'position',[560   200   600   700])
    suptitle(['Cluster=' num2str(clusters_to_consider(k)) ...
        ', n=' num2str(model.NA(clusters_to_consider(k))) ' flies'])
    colormap(jet)
    fig_count=fig_count+1;
    
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
    
    %colormap(jet)
    jet2=jet;
    jet2(size(jet2,1)/2,:)=[0 0 0];
    jet2(size(jet2,1)/2+1,:)=[0 0 0];
    colormap(jet2);
    
    % Figure 4b/left
    figure(fig_count)
    diff_state_normalized=NaN(length(bins),length(bins),model.Qdim);
    for i=1:model.Qdim
        diff_state_normalized(:,:,i) =during_state_normalized(:,:,i)-before_state_normalized(:,:,i);
        subplot(ceil(model.Qdim/3),3,i)
        imagesc(diff_state_normalized(:,:,i),[-1 1]);hold on;
        th = 0:pi/50:2*pi;
        R = length(bins)./2;
        xunit = R*1.5/3.2 * cos(th) + R;
        yunit = R*1.5/3.2 * sin(th) + R;
        plot(xunit, yunit,'w','LineWidth',1.5);
        xunit = R * cos(th) + R;
        yunit = R * sin(th) + R;
        plot(xunit, yunit,'w','LineWidth',1.5);
        hold off;
        
        axis off
    end
    set(gcf,'position',[560   200   600   700])
    suptitle(['Cluster=' num2str(clusters_to_consider(k)) ...
        ', n=' num2str(model.NA(clusters_to_consider(k))) ' flies'])
%     subplot(ceil(model.Qdim/3),3,i+1)
%     imagesc(diff_state_normalized(:,:,i),[-1 1]);hold on
    
    axis off
    colormap(jet2)
    colorbar
    %print('-dpdf',[folder_PDF 'Fig4b1'])
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
    hold off;
    fig_count = fig_count+1;
    
    before_norm{k} = before_state_normalized;
    during_norm{k} = during_state_normalized;
    diff_norm{k} = diff_state_normalized;
end

end
function [] = HLSDensitySpaceCircle(model,before_state_normalized,during_state_normalized,clusters_to_consider,fig_title)
%% plotting density of states as a function of space
if exist('fig_count','var')==0
    fig_count=90;
else
    fig_count=fig_count+1;
end

sz=size(before_state_normalized{1}(:,:,1));
x2=ceil(sz(1)/2); y2=ceil(sz(2)/2); % circular region parameters
[xgrid, ygrid] = meshgrid(1:sz(2), 1:sz(1));
x = xgrid - x2;    % offset the origin
y = ygrid - y2;
InsideMask = cell(1,11);
for i = 1:11
    r=(i-1)/10*(x2-1);
    % Make a logical image with the selected circular region set to 1, the rest
    % to zero
    InsideMask{i} = x.^2 + y.^2 <= r.^2;
end
InsideMask{1}(InsideMask{1}) = false;

for k=1:length(clusters_to_consider)
    plot_before = cell(model.Qdim,10);
    plot_during = cell(model.Qdim,10);
    for i = 1:model.Qdim
        img_before = before_state_normalized{k}(:,:,i);
        img_during = during_state_normalized{k}(:,:,i);
        for j=2:11
            r_o=(j-1)/10*(x2-1);
            r_i=(j-2)/10*(x2-1);
            plot_before{i,j-1} = sum(sum(double(img_before).*(InsideMask{j}-InsideMask{j-1})))./(pi*(r_o.^2-r_i.^2));
            plot_during{i,j-1} = sum(sum(double(img_during).*(InsideMask{j}-InsideMask{j-1})))./(pi*(r_o.^2-r_i.^2));
            
        end
    end
    
    r_val = [0:0.1:0.9,1.1];
    for i = 1:model.Qdim
        before_val = cell2mat(plot_before(i,:));
        during_val = cell2mat(plot_during(i,:));
        figure(fig_count);subplot(ceil(model.Qdim/3),3,i);stairs(r_val,[before_val,0]);
        hold on;stairs(-1*r_val,[before_val,0]);axis([-1,1,0,0.4]);
        line([1.5/3.2 1.5/3.2],[0 0.4]);line([-1.5/3.2 -1.5/3.2],[0 0.4]);hold off
        title(['HLS ' num2str(i) ' before'])
        xlabel('Normalized Distance')
        figure(fig_count+1);subplot(ceil(model.Qdim/3),3,i);stairs(r_val,[during_val,0])
        hold on;stairs(-1*r_val,[during_val,0]);axis([-1,1,0,0.4]);
        line([1.5/3.2 1.5/3.2],[0 0.4]);line([-1.5/3.2 -1.5/3.2],[0 0.4]);hold off
        title(['HLS ' num2str(i) ' during'])
        xlabel('Normalized Distance')
        
        figure(fig_count+2);subplot(ceil(model.Qdim/3),3,i);stairs(r_val,[during_val-before_val,0],'Color',[0 0 0])
        hold on;stairs(-1*r_val,[during_val-before_val,0],'Color',[0 0 0]);axis([-1,1,-0.15,0.15]);
        line([1.5/3.2 1.5/3.2],[-0.2 0.2],'color','k','linestyle','--');line([-1.5/3.2 -1.5/3.2],[-0.2 0.2],'color','k','linestyle','--');
        line([-1 1],[0 0],'color','k');hold off
        title(['HLS ' num2str(i) ' difference'])
        xlabel('Normalized Distance')
    end
    figure(fig_count);suptitle(['Cluster=' num2str(clusters_to_consider(k)) ...
        ', n=' num2str(model.NA(clusters_to_consider(k))) ' flies'])
    set(gcf,'position',[849 49 824 918])
    figure(fig_count+1);suptitle(['Cluster=' num2str(clusters_to_consider(k)) ...
        ', n=' num2str(model.NA(clusters_to_consider(k))) ' flies'])
    set(gcf,'position',[849 49 824 918])
    figure(fig_count+2);suptitle(['Cluster=' num2str(clusters_to_consider(k)) ...
        ', n=' num2str(model.NA(clusters_to_consider(k))) ' flies'])
    set(gcf,'position',[849 49 824 918])
    if ~isempty(fig_title)
        figure(fig_count);print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        figure(fig_count+1);print('-dpsc2',[fig_title '.ps'],'-loose','-append');
        figure(fig_count+2);print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
    fig_count = fig_count+3;
end


end

% Supplementary Figures
function [] = TransposeProperties(model,clusters_to_consider,ndx,ndxLL,fig_title)
%% plotting the basic model TPs (Figure 3)
if exist('fig_count','var')==0
    fig_count=1;
else
    fig_count=fig_count+1;
end
clims = [0 1];
for j=1:length(clusters_to_consider)
    figure(fig_count)
    subplot((ceil(model.Qdim/3)),4,[2 3 6 7]);
    TP = model.HHMMs{1,clusters_to_consider(j)}.Q.alpha;                    % TP is the higher level Transpose Property matrix
    TP_sorted=TP(ndx{j},ndx{j});                                            % sort TP based on ndx
    
    TP_tot = repmat(sum(TP_sorted,1),length(ndx{j}),1);
    TP_normalized = TP_sorted./TP_tot;
    
    imagesc(TP_normalized,clims);
    set(gca,'color',[1 1 1]);
    colorbar
    %sbpl = [1 5 9 13 14 15 16 12 8 4];
    sbpl = setdiff(1:1:(ceil(model.Qdim/3))*4,[2 3 6 7]);
    for jj=1:max(ndx{j})
        subplot((ceil(model.Qdim/3)),4,sbpl(jj))
        TP_LL = model.HHMMs{1,clusters_to_consider(j)}.A{1,ndx{j}(jj)}.alpha;
        TP_LL_sorted=TP_LL(squeeze(ndxLL{j,jj})',squeeze(ndxLL{j,jj})');
        if ~isempty(TP_LL_sorted)
            TP_LL_tot = repmat(sum(TP_LL_sorted,1),model.dim,1);
            TP_LL_normalized = TP_LL_sorted./TP_LL_tot;
            
            imagesc(TP_LL_normalized,clims)
            set(gca,'color',[1 1 1]);
            colorbar
        end
        title(['HLS #: ' num2str(jj)])
    end
    set(gcf, 'Position', [350   450   900   460])
    pause(0.5)
    suptitle(['Cluster ' num2str(clusters_to_consider(j)) ', ' num2str(model.NA(clusters_to_consider(j))) ' Flies'])
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
    fig_count=fig_count+1;
end
end


