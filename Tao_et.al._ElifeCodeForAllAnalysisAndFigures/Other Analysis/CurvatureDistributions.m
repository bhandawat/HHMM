function [] = CurvatureDistributions(model,HHMM_by_State_sorted,fig_title)
nHLS = model.Qdim;
tmp = squeeze(HHMM_by_State_sorted);

allCurv = cell(size(tmp,2),1);allCurvNorm = cell(size(tmp,2),1);
for j = 1:size(tmp,2)
    tmp2 = squeeze(tmp(:,j,:));
    tmp2 = tmp2(~cellfun('isempty',tmp2));
    kin = cellfun(@(x) sum(x,2),tmp2,'UniformOutput', false);
    len = cell2mat(cellfun(@(x) size(x,2),tmp2,'UniformOutput', false));
    curv = cellfun(@(x) x(16),kin);
    
    curv(len<20) = [];
    kin(len<20) = [];
    len(len<20) = [];
    allCurv{j} = curv;
    allCurvNorm{j} = curv./len.*30;
end

cc=varycolor(nHLS+1);
cc(end,:) = [1 1 1];
x_values = linspace(-20,20,100);
y = zeros(nHLS,100);
for i = 1:nHLS
    [pdca] = fitdist(allCurvNorm{i},'kernel');
    y(i,:) = pdf(pdca,x_values);
end
figure;hold on
for i = 8:nHLS
    plot(x_values,y(i,:),'Color',cc(i,:))
end
plot([0 0],[0 0.14],'k');
text(-18,0.13,['HL State 8 med: ' num2str(median(allCurvNorm{8})) ' rad/s'])
text(-18,0.12,['HL State 9 med: ' num2str(median(allCurvNorm{9})) ' rad/s'])
text(-18,0.11,['HL State 10 med: ' num2str(median(allCurvNorm{10})) ' rad/s'])

text(-18,0.10,['totTurn 8: ' num2str(sum(allCurv{8})) ' rad'])
text(-18,0.09,['totTurn 9: ' num2str(sum(allCurv{9})) ' rad'])
text(-18,0.08,['totTurn 10: ' num2str(sum(allCurv{10})) ' rad'])

ylabel('PDF');xlabel('Track curvature (radians/sec)')
legend({'HL State 8','HL State 9','HL State 10','zero line'})
title('Tracks >= 10 frames')

if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

end