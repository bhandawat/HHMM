function [] = simplexXMeans(model,options,fig_title)
if isempty(options.empKMFile)
    disp('need kmeans GLM file')
end
load(options.empKMFile,'kMeans')

kMeansGlob = cell(1,length(kMeans));
for c = 1:length(kMeans)
    nFlys = size(kMeans{c}.Seed,2);
    features = num2cell(drchrnd(ones(1,model.Qdim), nFlys),2);
    distance = 'euclidean';
    maxCluster = 0;
    
    for ndx = 1:nFlys
        kMeansTemp = ml_kmeans_aicbic(features, [], distance, ndx, 1);
        diff_bic=-diff(kMeansTemp.bic);
        df = length(features)-1;
        t = sqrt(log(df+1)+diff_bic);
        t_val = real(t);
        P_val = 2*tcdf(-abs(t_val),df);
        perCond = (kMeansTemp.bic-min(kMeansTemp.bic))./min(kMeansTemp.bic)<0.05;
        optCluster = find(P_val<0.05 & perCond(1:end-1),1,'last')+1;
        
        kMeansTemp.optCluster = optCluster;
        kMeansTemp.P_val = [NaN P_val];
        kMeansTemp.t_val = [NaN t_val];
        
        
        if ndx == 1
            kMeansGlob{c} = kMeansTemp;
        end
        
        if ~isempty(optCluster)
            if optCluster>maxCluster
                maxCluster = optCluster;
                kMeansGlob{c} = kMeansTemp;
            end
        end
    end
end

for c = 1:length(kMeans)
    figure(c);
    k = 1;
    for state = 1:2
        for loc = 1:2
            optC = kMeans{c}.maxCluster.data{loc,state}.optCluster;
            BIC = kMeans{c}.maxCluster.data{loc,state}.bic;
            AIC = kMeans{c}.maxCluster.data{loc,state}.aic;
            subplot(2,2,k);
            plot(BIC,'Linewidth',2);hold on
            plot(kMeansGlob{c}.bic,'Linewidth',2);
            plot(optC,BIC(optC),'ro');
            k = k+1;
            legend({'Empirical BIC','Rand BIC'},'Location','northWest')
            title([kMeans{c}.scenario{loc,state} ', Opt C: ' num2str(optC)])
            xlabel('k clusters');
            ylabel('BIC score')
            hold off
        end
    end
    set(gcf,'position',[849 49 824 918])
    if ~isempty(fig_title)
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
end

figure;
for c = 1:length(kMeansGlob)
    if length(kMeansGlob)>2
        subplot(ceil(length(kMeansGlob)./2),2,c)
    else
        subplot(2,1,c)
    end
    plot(kMeansGlob{c}.bic,'Linewidth',2);hold on
    plot(kMeansGlob{c}.aic,'Linewidth',2);
    legend({'bic','aic'},'Location','northWest')
end
if ~isempty(fig_title)
    print('-dpsc2',[fig_title '.ps'],'-loose','-append');
end

end