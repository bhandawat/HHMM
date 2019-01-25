function [] = forwardTracks(longTrack,maxTrack)
ColorSet=varycolor(100);
states_of_interest = find(maxTrack>0);
kk = 1;
figure;set(gcf,'position',[9 49 824 918])
for i=1:length(states_of_interest)
    state = states_of_interest(i);
    allTrack = longTrack(state,:);
    allTrack = cat(2, allTrack{:});
    
    h1=subplot(4,3,kk);
    set(h1, 'ColorOrder', ColorSet);
    hold all;
    
    if ~isempty(allTrack)
        allTrack = allTrack(~cellfun('isempty',allTrack));
        for k=1:length(allTrack)
            x = allTrack{k}(1,:);
            y = allTrack{k}(2,:);
            x = x-x(1);                                                     % make everything relative to the first point of trajectory (x)
            y = y-y(1);                                                     % (y)
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
        title(['State ' num2str(state)])
        kk = kk+1;
    end
    
    if kk >25
        figure;set(gcf,'position',[9 49 824 918])
        kk = 1;
    end
end

end