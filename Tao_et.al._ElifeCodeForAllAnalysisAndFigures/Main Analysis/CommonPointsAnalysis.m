function [] = CommonPointsAnalysis(commTracksAll,data_time,fig_title)
close all
for c = 1:size(commTracksAll,1)
    if ~isempty(commTracksAll{c,1})
        commTracks = commTracksAll(c,:);
        for i = 1:2
            for j = 1:size(commTracks{i},1)
                d = data_time{commTracks{i}(j,3)};
                track{j,i} = d(:,max(commTracks{i}(j,1),11):commTracks{i}(j,2));
                if commTracks{i}(j,1)<1
                    track2{j,i} = [];
                else
                    track2{j,i} = d(:,max(commTracks{i}(j,1),11)-10:commTracks{i}(j,2));
                end
            end
        end
        forwardTrack(track2)
        VParVperp(track,commTracks,fig_title)
        close all
    end
end

end

function [] = forwardTrack(track)

figure(3)
ColorSet=varycolor(100);

for i = 1:size(track,2)
    h1=subplot(3,2,i);
    hold on
    set(h1, 'ColorOrder', ColorSet);
    
    currTrack = track(:,i);
    currTrack = currTrack(~cellfun('isempty',currTrack));
    
    for j = 1:size(currTrack,1)
        x = currTrack{j,1}(1,:);
        y = currTrack{j,1}(2,:);
        
        if length(x)>11
            x = x-x(1);
            y = y-y(1);
            xy  = x+1i*y;                                                   % change to polar coordinates for easier calculation
            vel_x=diff(x); vel_y=diff(y);
            vel=vel_x+vel_y*1i;
            theta=angle(sum(vel(1:10)));                                     % assume that the first 10 steps define the curvature (i.e 1 and 8 are both in the direct forward direction)
            xy2=xy.*exp(-theta*1i+(pi/2)*1i);                               % shift by 90 degrees to point up instead of to the right
            newY = imag(xy2);
            if(length(newY(newY>0))>length(newY(newY<0))) && (max(abs(sqrt(vel_x.^2+vel_y.^2)))<0.7)
                plot(real(xy2),imag(xy2));                                  % only plot tracks where the data doesn't jump around too much
            end
        end
    end
    xlim([-0.1 0.1])
    ylim([-0.05 0.25])
    title(['ForwardTracks, HLS: ' num2str(i) ', n= ' num2str(j)])
end

for i = 1:size(track,2)
    subplot(3,2,i+2);
    currTrack = track(:,i);
    currTrack = currTrack(~cellfun('isempty',currTrack));
    for j = 1:size(currTrack,1)
        vPar = currTrack{j,1}(15,:);
        vPerp = currTrack{j,1}(16,:);
        plot(vPar(11:end),vPerp(11:end),'Color',[0.7 0.7 0.7])
        hold on
    end
    title(['PhasePlane, HLS: ' num2str(i) ', n= ' num2str(j)])
    xlabel('vPar');ylabel('vPerp')
    xlim([-0.1 0.7]);ylim([-0.25 0.25])
    hold off
end


end

function [] = VParVperp(track,commTracks,fig_title)

vParAll = cell(2,1);vPerpAll = cell(2,1);
frameWin = 5;
for i = 1:size(track,2)
    
    currTrack = track(:,i);
    currTrack = currTrack(~cellfun('isempty',currTrack));
    ndxFromBeg = commTracks{i}(:,4)-commTracks{i}(:,1);
    ndxFromEnd = commTracks{i}(:,2)-commTracks{i}(:,4);
    kk = 1;
    for j = 1:size(ndxFromBeg,1)
        if (ndxFromBeg(j) >=frameWin) && (ndxFromEnd(j) >=frameWin)
            vPar = track{j,i}(15,:);
            vPerp = track{j,i}(16,:);
            if length(vPar)>10
                vParAll{i}(kk,:) = vPar(ndxFromBeg(j)-frameWin+1:ndxFromBeg(j)+frameWin+1);
                vPerpAll{i}(kk,:) = vPerp(ndxFromBeg(j)-frameWin+1:ndxFromBeg(j)+frameWin+1);
            end
            kk = kk+1;
        end
    end
end

figure(1)
for i = 1:size(track,2)
    subplot(2,2,i)
    for j = 1:size(vParAll{i},1)
        plot(vParAll{i}(j,:),'Color',[0.7 0.7 0.7])
        hold on
    end
    plot(median(vParAll{i}),'k','LineWidth',1)
    title(['vPar, HLS: ' num2str(i) ', n= ' num2str(j)])
    xlabel('frames');ylabel('vPar')
    ylim([-0.1 0.7])
    hold off
end

for i = 1:size(track,2)
    subplot(2,2,i+2)
    for j = 1:size(vPerpAll{i},1)
        plot(vPerpAll{i}(j,:),'Color',[0.7 0.7 0.7])
        hold on
    end
    plot(median(vPerpAll{i}),'k','LineWidth',2)
    title(['vPerp, HLS: ' num2str(i) ', n= ' num2str(j)])
    xlabel('frames');ylabel('vPerp')
    ylim([-0.25 0.25])
    hold off
end
set(gcf,'position',[9 49 824 918])
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end


figure(3)
for i = 1:size(track,2)
    subplot(3,2,i+4)
    for j = 1:size(vPerpAll{i},1)
        plot(vParAll{i}(j,:),vPerpAll{i}(j,:),'Color',[0.7 0.7 0.7])
        hold on
    end
    %plot(median(vParAll{i}),median(vPerpAll{i}),'k','LineWidth',2)
    title(['PhasePlane, HLS: ' num2str(i) ', n= ' num2str(j)])
    xlabel('vPar');ylabel('vPerp')
    xlim([-0.1 0.7]);ylim([-0.25 0.25])
    hold off
end
set(gcf,'position',[9 49 824 918])
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

% figure(4)
% plot(vParAll{1}(:,frameWin+1),vPerpAll{1}(:,frameWin+1),'.');hold on
% plot(vParAll{2}(:,frameWin+1),vPerpAll{2}(:,frameWin+1),'.');
% xlim([-0.2 6]);ylim([-2 2])

end