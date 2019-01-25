function [] = TemporalAnalysis(data_time,figName)

firstEntry = reshape(zeros(size(data_time)),[],1);
nFly = size(data_time,2);
for i = 1:nFly
    rPos = sqrt(data_time{i}(1,:).^2+data_time{i}(2,:).^2);
    Out = rPos>(1.9/3.2);
    During = [false(1,floor(length(rPos)/2)), true(1,ceil(length(rPos)-length(rPos)/2))];

    firstEntry(i) = find(During & ~Out,1);% odor on and fly inside
end

% radial probability
figNum = 10;
tBin = [0 30 150 300 5400];
Bins = cell(1,length(tBin)-1);
dx = 0.05;edges = 0:dx:1;border = 1.5/3.2;
c=gray(length(tBin));
leg = cell(1,length(tBin)-1);
figure(figNum);set(gcf,'position',[855 49 824 918])
for i = 1:length(tBin)-1
    for fly = 1:nFly
        rPos = sqrt(data_time{fly}(1,:).^2+data_time{fly}(2,:).^2);
        Bins{i} = [Bins{i} rPos(firstEntry(fly)+tBin(i):min(firstEntry(fly)+tBin(i+1),length(rPos)))];
    end
    [N,~] = histcounts(Bins{i},edges);
    N = N+eps;
    Prob.y{i} = N./sum(N);
    Prob.x{i} = edges(2:end)-dx./2;

    subplot(2,1,1);
    plot(Prob.x{i},Prob.y{i},'Color',c(i,:),'LineWidth',2);hold on;
    leg{i} = [num2str(tBin(i)./30) ' to ' num2str(tBin(i+1)./30) ' s'];
end
plot([border,border],[0 1],'r--');
xlabel('Radial Distance');ylabel('Probability');
title('Radial Occupancy')
legend(leg)
xlim([0 1]);ylim([0 1])


% radial probability
tBin = [0 300 600 900 1200 1500 5400];
Bins = cell(1,length(tBin)-1);
c=gray(length(tBin));
leg = cell(1,length(tBin)-1);
for i = 1:length(tBin)-1
    for fly = 1:nFly
        rPos = sqrt(data_time{fly}(1,:).^2+data_time{fly}(2,:).^2);
        Bins{i} = [Bins{i} rPos(firstEntry(fly)+tBin(i):min(firstEntry(fly)+tBin(i+1),length(rPos)))];
    end
    [N,~] = histcounts(Bins{i},edges);
    N = N+eps;
    Prob.y{i} = N./sum(N);
    Prob.x{i} = edges(2:end)-dx./2;

    subplot(2,1,2);
    plot(Prob.x{i},Prob.y{i},'Color',c(i,:),'LineWidth',2);hold on
    leg{i} = [num2str(tBin(i)./30) ' to ' num2str(tBin(i+1)./30) ' s'];
end
plot([border,border],[0 1],'r--');
xlabel('Radial Distance');ylabel('Probability');
title('Radial Occupancy')
legend(leg)
xlim([0 1]);ylim([0 1])

if ~isempty(figName)
    for i = 1:figNum
        figure(i)
        print('-painters','-dpsc2',[figName '.ps'],'-loose','-append');
    end
end

end

