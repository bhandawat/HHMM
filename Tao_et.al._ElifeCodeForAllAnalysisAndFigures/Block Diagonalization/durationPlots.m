function [] = durationPlots(dur_HHMM_sorted,dur_HMM_States,dur_HMM_sorted,statesOfInterest,HHMM_lowlevel_sorted,fig_title)

figure(10);
temp = cell2mat(dur_HMM_States');
n = histcounts(temp,1:max(temp)+1);
n = n.*(1:1:(max(temp)));
ccdf = 1-cumsum(n)./sum(n);
subplot(3,1,1);hold on;plot((1:1:max(temp))./30,ccdf,'b')
xlabel('Time (s)')
ylabel('% data')
xlim([0 3]);ylim([0 1])
text(1.5,0.6, ['Total # frames: ' num2str(sum(temp))],'Color','b')
text(1.5,0.8, ['Total # of tracks: ' num2str(length(temp))],'Color','b')

figure(10);
for i = 2:length(statesOfInterest)
    temp = cell2mat(dur_HMM_sorted(statesOfInterest{i})');
    n = histcounts(temp,1:max(temp)+1);
    n = n.*(1:1:(max(temp)));
    ccdf = 1-cumsum(n)./sum(n);
    subplot(3,1,i);hold on;plot((1:1:max(temp))./30,ccdf,'b')
    xlabel('Time (s)')
    ylabel('% data')
    xlim([0 3]);ylim([0 1])
    title(['States: ' num2str(statesOfInterest{i})])
    text(1.5,0.6, ['Total # frames: ' num2str(sum(temp))],'Color','b')
    text(1.5,0.8, ['Total # tracks: ' num2str(length(temp))],'Color','b')
end

statesOfInterest = {[1:10],1,10};
for i = 1:length(statesOfInterest)
    temp = cell2mat(dur_HHMM_sorted(statesOfInterest{i})');
    n = histcounts(temp,1:max(temp)+1);
    n = n.*(1:1:(max(temp)));
    ccdf = 1-cumsum(n)./sum(n);
    subplot(3,1,i);hold on;plot((1:1:max(temp))./30,ccdf,'r')
    legend({'HMM','HHMM'})
    xlabel('Time (s)')
    ylabel('% data')
    xlim([0 3]);ylim([0 1])
    title(['States: ' num2str(statesOfInterest{i})])
    text(1.5,0.5, ['Total # frames: ' num2str(sum(temp))],'Color','r')
    text(1.5,0.7, ['Total # tracks: ' num2str(length(temp))],'Color','r')
end

subplot(3,1,1)
title('States: all')
subplot(3,1,2)
title('States: 3,1')
subplot(3,1,3)
title('States: 7,10')
set(gcf,'position',[849 49 824 918])
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

plottingLLSHLSDistributions(dur_HHMM_sorted,HHMM_lowlevel_sorted,fig_title);

end


function [] = plottingLLSHLSDistributions(dur_HHMM_sorted,HHMM_lowlevel_sorted,fig_title)

nHLS = size(HHMM_lowlevel_sorted,3);

LLS_dur = cell(nHLS,3000);
LLS_nTrans = cell(nHLS,3000);
for i = 1:size(HHMM_lowlevel_sorted,1)
    for j = 1:nHLS
        currHL = squeeze(HHMM_lowlevel_sorted(i,:,j,:));
        currHL = currHL(cellfun(@numel,currHL)>0);
        for k = 1:length(currHL)
            [~,currLLS] = max(currHL{k}(:,11:end));
            ii = [0, diff(currLLS(:)')==0,0];
            durTmp = diff(find([1,diff(currLLS),1]));
            nTransTmp = length(durTmp)-1;
            
            LLS_dur{j,length(currLLS)} = [LLS_dur{j,length(currLLS)} durTmp];
            LLS_nTrans{j,length(currLLS)} = [LLS_nTrans{j,length(currLLS)} nTransTmp];
        end
        
    end
end

bins = [1:1:20 30:10:60];
avgNTrans = zeros(nHLS,max(bins));avgNTrans2 = zeros(nHLS,1);
avgDur = zeros(nHLS,max(bins));
for i = 1:nHLS
    for j = 1:length(bins)-1
        xVals = bins(j):bins(j+1)-1;
        avgNTrans(i,xVals) = mean(cell2mat(LLS_nTrans(i,xVals)));
        avgDur(i,xVals) = mean(cell2mat(LLS_dur(i,xVals)));
    end
    avgNTrans2(i) = mean(cell2mat(LLS_nTrans(i,:)));
end

for j = 1:length(bins)-1
    xVals = bins(j):bins(j+1)-1;
    avgNTrans(nHLS+1,xVals) = mean(cell2mat(reshape(LLS_nTrans(:,xVals),1,[])));
    avgDur(nHLS+1,xVals) = mean(cell2mat(reshape(LLS_dur(:,xVals),1,[])));
    avgNTrans2(nHLS+1) = mean(cell2mat(reshape(LLS_nTrans,1,[])));
end


% plot ccdf of tracks for LLS and HLS
figure;set(gcf,'position',[2 42 838 924])
dur_HHMM_sorted{11} = cell2mat(dur_HHMM_sorted(1:10)');
LLS_dur{11,2} = cell2mat(reshape(LLS_dur(1:nHLS,:),1,[]));
for i = 1:nHLS+1
    temp = cell2mat(LLS_dur(i,:));
    n = histcounts(temp,1:max(temp)+1);
    ccdf = 1-cumsum(n)./sum(n);
    
    tempHLS = dur_HHMM_sorted{i};
    durHLS = histcounts(tempHLS,1:max(tempHLS)+1);
    ccdfHLS = 1-cumsum(durHLS)./sum(durHLS);
    
    subplot(4,3,i);yyaxis left
    plot((0:1:min(length(ccdf),45))./30,[1 ccdf(1:1:min(length(ccdf),45))])
    hold on;plot((0:1:min(length(ccdfHLS),45))./30,[1 ccdfHLS(1:1:min(length(ccdfHLS),45))])
    xlabel('Time (s)')
    ylabel('% data')
    ylim([0 1]);
    
    yyaxis right
    stairs((1:45)./30,avgNTrans(i,1:45));
    ylabel('Mean # of transitions')
    xlim([0 1.5]);ylim([0 6])
    
    text(0.2, 5, ['Avg # Trans = ' num2str(avgNTrans2(i))])
    title(['HLS: ' num2str(i)])
end
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end


% plot average number of LLS segments/HLS track
figure;set(gcf,'position',[2 42 838 924])
for i = 1:nHLS+1
    subplot(4,3,i);
    stairs((1:45)./30,avgDur(i,1:45)./30,'r')
    ylim([0 0.6]);xlim([0 1.5])
    xlabel('Time (s)');ylabel('Mean duration of LLS')
end
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

cc=varycolor(nHLS+1);
cc(end,:) = 0;cc = [cc;zeros(1,3)];
figure;set(gcf,'position',[2 42 838 924]);hold on;
dur_HHMM_sorted{11} = cell2mat(dur_HHMM_sorted(1:10)');
dur_HHMM_sorted{12} = cell2mat(dur_HHMM_sorted(3:10)');
leg = cell(1,11);
for i = 1:nHLS+2
    temp = dur_HHMM_sorted{i};
    n = histcounts(temp,1:max(temp)+1);
    n = n.*(1:1:(max(temp)));
    ccdf = 1-cumsum(n)./sum(n);
    if i == 12
        plot((1:min(length(ccdf),90))./30,ccdf(1:min(length(ccdf),90)),'--','Color',cc(i,:))
    else
        plot((1:min(length(ccdf),90))./30,ccdf(1:min(length(ccdf),90)),'Color',cc(i,:))
    end
    leg{i} = ['HLS ' num2str(i)];
end
leg{nHLS+1} = 'All HLS';
leg{nHLS+2} = 'HLS 3 to 10';
xlabel('Time (s)')
ylabel('% data')
xlim([0 3]);ylim([0 1])
title('Duration Distribution (All points)')
legend(leg)
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end




end

