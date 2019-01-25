function [] = SpaceTimeAnalysis(data_time,model,likely_high_state_sorted,c2cons,figName)
nHLS = model.Qdim;
nPts = ceil(size(data_time{1},2)./2);
RThresh = 1.9/3.2;BinThresh = ceil(RThresh.*100);

for clust = 1:length(c2cons)
    nFlys = length(c2cons{clust});
    if nFlys>0
        close all
        sDI = zeros(nPts,100,nHLS);sBI = zeros(nPts,100,nHLS);
        sDO = zeros(nPts,100,nHLS);sBO = zeros(nPts,100,nHLS);
        for i = 1:nFlys
            c2cons{clust}(i);
            r = sqrt(data_time{1,i}(1,:).^2+data_time{1,i}(2,:).^2);
            half = floor(length(r)./2);
            idxDI = [false(1,half) (r(half+1:end)<RThresh)];
            idxBI = [(r(1:half)<RThresh) false(1,half+1)];
            idxDO = [false(1,half) (r(half+1:end)>=RThresh)];
            idxBO = [(r(1:half)>=RThresh) false(1,half+1)];
            
            [sDI] = calcSpaceTime(sDI,r(idxDI),i,idxDI(idxDI),likely_high_state_sorted(:,idxDI));
            [sDO] = calcSpaceTime(sDO,r(idxDO),i,idxDO(idxDO),likely_high_state_sorted(:,idxDO));
            [sBI] = calcSpaceTime(sBI,r(idxBI),i,idxBI(idxBI),likely_high_state_sorted(:,idxBI));
            [sBO] = calcSpaceTime(sBO,r(idxBO),i,idxBO(idxBO),likely_high_state_sorted(:,idxBO));
        end
        
        sI = squeeze(sum(sBI(:,1:BinThresh,:),2));
        sI(:,:,2) = squeeze(sum(sDI(:,1:BinThresh,:),2));
        sO = squeeze(sum(sBO(:,BinThresh:end,:),2));
        sO(:,:,2) = squeeze(sum(sDO(:,BinThresh:end,:),2));
        
        % downsample Time
        t = 60;
        sI2 = zeros(size(sI,1)./t,size(sI,2),size(sI,3));
        sO2 = zeros(size(sO,1)./t,size(sO,2),size(sO,3));
        for i = 1:t
            sI2 = sI2+sI(i:t:end,:,:);
            sO2 = sO2+sO(i:t:end,:,:);
        end
        
        sINorm = sI2./sum(sI2,2);sINorm(isnan(sINorm)) = 0;
        sONorm = sO2./sum(sO2,2);sONorm(isnan(sONorm)) = 0;
        
        % plot inside (all HL on same plot During)
        c = varycolor(nHLS+1);leg = cell(1,nHLS);
        figure;
        subplot(2,1,1);set(gcf,'Position',[2 42 838 924])
        hold on
        for HL = 1:nHLS
            plot((1:t:size(sI,1))./30,sINorm(:,HL,2),'Color',c(HL,:));
            leg{HL} = ['HLS ' num2str(HL)];
        end
        ylim([0 0.6]);xlim([0 60])
        legend(leg,'Location','Best')
        xlabel('time since entry (s)')
        ylabel('Probability')
        suptitle('P(HL state| time) Inside')
        
        % plot inside (Grouped HL on same plot During)
        g = [{1},{2},{3},{4:7},{8:10}];
        c = varycolor(length(g));leg = cell(1,length(g));
        subplot(2,1,2);set(gcf,'Position',[2 42 838 924])
        hold on
        for HL = 1:length(g)
            plot((1:t:size(sI,1))./30,sum(sINorm(:,g{HL},2),2)./length(g{HL}),'Color',c(HL,:));
            leg{HL} = ['HLS ' num2str(g{HL})];
        end
        ylim([0 0.6]);xlim([0 60])
        legend(leg,'Location','Best')
        xlabel('time since entry (s)')
        ylabel('Probability')
        suptitle('P(HL state| time) Inside')
        if ~isempty(figName)
            print('-painters','-dpsc2',[figName '.ps'],'-loose','-append');
        end
    end
end

end

function [s] = calcSpaceTime(s,r,fly,idx,likely_high_state_sorted)
ii = [0, diff(idx(:)')==0,0];
i1 = strfind(ii,[0 1]);
i2 = strfind(ii,[1 0]);

for j = 1:length(i1)
    if idx(i1(j))
        t = (i1(j):i2(j));
        d = ceil(r(t).*100);
        HLS = likely_high_state_sorted(fly,t);
        t2 = t-t(1)+1;
        
        t2(HLS == 0) = [];
        d(HLS == 0) = [];d(d>100) = 100;
        HLS(HLS == 0) = [];
        
        for k = 1:length(t2)
            s(t2(k),d(k),HLS(k)) = s(t2(k),d(k),HLS(k))+1;
        end
    end
end
end








