function [] = bootstrapping(model,B,HLSDistRawOn,HLSDistRawOff,bad_Fly_List)
close all
goodfly = cell(1,2);
BNew = cell(2,3);BNew2 = cell(2,2);
scenarios = {'during in','during out','before in','before out','diff in','diff out'};
for l = 1:2
    goodfly{l} = 1:34;
    goodfly{l}(bad_Fly_List) = [];
    B1 = cat(2,HLSDistRawOn{goodfly{l},l});
    B2 = cat(2,HLSDistRawOff{goodfly{l},l});
    BNew2(l,:) = {sum(B1,2),sum(B2,2)};
    
    B1 = sum(B1,2)./sum(sum(B1,2),1);
    B2 = sum(B2,2)./sum(sum(B2,2),1);
    
    B3 = B1-B2;
    BNew(l,:) = {B1,B2,B3};
    
end

nIt = 10000;
nSamp = 34;
CI = cell(size(BNew2));
k = 1;
figure(1);set(gcf,'position',[849 49 824 918])
for state = 1:2
    for loc = 1:2
        xDist = zeros(nIt,model.Qdim);
        
        for HLS = 1:model.Qdim
            p = B{loc,state}(HLS,:);
            x = reshape(discretesample2(ones(1,length(p))./length(p), nSamp.*nIt),nIt,nSamp);
            %[CI{loc,state}(:,HLS),avg] = calculateCI(mean(p(x),2));
            CI{loc,state}(:,HLS) = bootci(nIt,{@mean, p});
        end
        
        meanBootStrap = (CI{loc,state}(2,:)+CI{loc,state}(1,:))./2;
        errBootStrap = (CI{loc,state}(2,:)-CI{loc,state}(1,:))./2;
        
        figure(1)
        subplot(4,2,k)
        bar(meanBootStrap);hold on;
        errorbar(1:1:model.Qdim,meanBootStrap,errBootStrap,'.')
        xlim([0 11]);title(scenarios{k})
        if state ==3
            ylim([-0.2 0.2])
        else
            ylim([0 0.5])
        end
        k = k+1;
        
        
    end
end

p_val = zeros(2,model.Qdim);
for loc = 1:2
    for HLS = 1:model.Qdim
        p = [B{loc,1}(HLS,:),B{loc,2}(HLS,:)];
        [~,bootNdx] = bootstrp(nIt,@mean,p);
        bootsam = p(bootNdx);
        xDur = bootsam(1:nSamp,:);
        xBef = bootsam(nSamp+1:end,:);
        
        t_obs = calculateT(p(1:nSamp)',p(nSamp+1:end)');
        if t_obs<0
            t_obs=-t_obs;
            t_prime = calculateT(xBef,xDur);
        else
            t_prime = calculateT(xDur,xBef);
        end
        p_val(loc,HLS) = sum(t_prime>t_obs)./nIt;
    end
end


p_val2 = zeros(2,model.Qdim);
for loc = 1:2
    for HLS = 1:model.Qdim
        p = [B{loc,1}(HLS,:),B{loc,2}(HLS,:)];
        p2 = p;
        p2(1:nSamp) = p(1:nSamp)-mean(p(1:nSamp))+mean(p);
        p2(nSamp+1:end) = p(nSamp+1:end)-mean(p(nSamp+1:end))+mean(p);
        
        [~,bootNdx1] = bootstrp(nIt,@mean,p2(1:nSamp));
        [~,bootNdx2] = bootstrp(nIt,@mean,p2(nSamp+1:end));
        
        bootsam = p2([bootNdx1;bootNdx2]);
        xDur = p2(bootNdx1);
        xBef = p2(bootNdx1+nSamp);
        
        t_obs = calculateT_UnequalVar(p(1:nSamp)',p(nSamp+1:end)');
        if t_obs<0
            t_obs=-t_obs;
            t_prime = calculateT_UnequalVar(xBef,xDur);
        else
            t_prime = calculateT_UnequalVar(xDur,xBef);
        end
        p_val2(loc,HLS) = sum(t_prime>t_obs)./nIt;
    end
end



for loc = 1:2
    subplot(4,2,k)
    bar(mean(B{loc,3},2));hold on;
    for HLS = 1:model.Qdim
        text(HLS-0.5,0.12+0.02*mod(HLS,3),num2str(round(p_val(loc,HLS),3)),'FontSize',8);
    end
    xlim([0 11]);ylim([-0.2 0.2]);
    title('Difference')
    
    subplot(4,2,k+2)
    bar(mean(B{loc,3},2)./mean(B{loc,2},2))
    xlim([0 11]);ylim([-1 1]);
    title('Fractional Change')
    k = k+1;
end

end

function [t_val] = calculateT_Basic(x1,x2)
t_val = (mean(x1)-mean(x2));
end

function [t_val] = calculateT(x1,x2)
n = size(x1,1);m = size(x2,1);
sig = sqrt((sum((x1-mean(x1)).^2)+...
    sum((x2-mean(x2)).^2))./(n+m-2));
normStat = sig.*sqrt((1/n)+(1/m));

t_val = (mean(x1)-mean(x2))./normStat;
end

function [t_val] = calculateT_UnequalVar(x1,x2)
n = size(x1,1);m = size(x2,1);
sig_1 = sum((x1-mean(x1)).^2)./(n-1);
sig_2 = sum((x2-mean(x2)).^2)./(m-1);
normStat = sqrt((sig_1/n)+(sig_2/m));

t_val = (mean(x1)-mean(x2))./normStat;
end










