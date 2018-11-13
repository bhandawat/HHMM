clear all
close all

load WT_ACV0_May2017.mat
load data2
data=data2; clear data2;
load('wild_type_final_model.mat')
model.update(data,1,1);
dc=0:0.01:1;
[m,idx]=max(model.p);
cc=jet(model.Qdim);
plotidx=0;
opts=statset('glmfit');
opts.MaxIter=4000;

for i=1:model.NC
    XX{i}=[];
    XX2{i}=[];
    Y{i}=[];
    rawX{i}=[];
    rawY{i}=[];
    rawX2{i}=[];
    rawY2{i}=[];
    X{i}=[];
    Y{i}=[];
    X2{i}=[];
    Y2{i}=[];
    idx2=find(idx==i);  % idx of datapoints assigned to cluster i
    np=length(idx2);
    if(np>0)
        plotidx=plotidx+1;
        
        px=ceil(sqrt(np));
        py=ceil(np/px);
        for j=1:np
            d1=data{idx2(j)}(:,1:end-1);
            o1=d1(12,1:end-1);
            
            TP1=model.HHMMs{i}.getQxi(data(idx2(j)));
            TP1=TP1{1};
            TP1=reshape(TP1,model.Qdim^2,size(TP1,3));
            TPidx{i} = reshape(model.HHMMs{i}.Q.alpha>1,model.Qdim^2,1);
            
            TP1=TP1(TPidx{i},:);
            
            inside = d1(11,1:end-1);
            p1=model.HHMMs{i}.p{idx2(j)};
            p1=p1(:,1:end-1);
            clear p2
            for k=1:model.Qdim
                p2(k,:)=sum(p1(model.HHMMs{i}.Aidx{k},:),1);
            end
            [m2,idx3]=max(p2);
            figure(2*plotidx-1);        
            subplot(px,py,j), scatter(d1(1,:),d1(2,:),3*ones(size(d1(1,:))),cc(idx3,:))
            subplot(px,py,j), title(['Cluster ',num2str(i), 'with P = ',num2str(model.p(i,idx2(j)))])
             figure(2*plotidx)
%             idxoff=find(o1==0);
             idxoff=find(o1==0 & inside);
             idxon=find(o1==1 & inside);
%             idxon=idxon(1:min(1000,length(idxon)));
%             temp=find(diff(idxon)>1);
%             idxon=idxon(temp(1)+1:end);
             
             
%             idxoff=idxoff(idxoff>min(idxon));
             len=min(length(idxon),length(idxoff));
             idxoff=idxoff(1:len);             
             idxon=idxon(1:len);             
%             idxoff=idxoff(idxoff<5400); 
%             idxoff=idxoff(end-length(idxon)+1:end);
             pon=p2(1:end-1,idxon);
             poff=p2(1:end-1,idxoff);
                          
             TPon=TP1(1:end-1,idxon);
             TPoff=TP1(1:end-1,idxoff);

             subplot(np,3,3*j-2), bar(mean(p2(:,idxoff)')), title('odor off')
             subplot(np,3,3*j-1), bar(mean(p2(:,idxon)')), title('odor on')
             subplot(np,3,3*j),   bar(mean(p2(:,idxon)')-mean(p2(:,idxoff)')), title('difference')
             
             d2on=d1(7:8,idxon)';
             d2on(:,2)=abs(d2on(:,2));
             d2off=d1(7:8,idxoff)';
             d2off(:,2)=abs(d2off(:,2));
             
             pon=[pon;d2on'];  %THESE ARE THE APPENDING LINES
             poff=[poff;d2off'];
             
             
             bs=5;
             clear pon2 TPon2 d22on TPoff2 poff2 d22off
             for k=1:floor(size(pon,2)/bs)
                 pon2(:,k)=mean(pon(:,(k-1)*bs+1:k*bs)')';
                 TPon2(:,k)=mean(TPon(:,(k-1)*bs+1:k*bs)')';                 
                 d22on(k,:)=mean(d2on((k-1)*bs+1:k*bs,:));
             end
             for k=1:floor(size(poff,2)/bs)
                 poff2(:,k)=mean(poff(:,(k-1)*bs+1:k*bs)')';
                 d22off(k,:)=mean(d2off((k-1)*bs+1:k*bs,:));
                 TPoff2(:,k)=mean(TPoff(:,(k-1)*bs+1:k*bs)')';
             end

             type{i}.fly{j}.rawX=[d2off;d2on;];
             type{i}.fly{j}.rawY=[zeros(length(idxoff),1);ones(length(idxon),1)];
             
             type{i}.fly{j}.rawX2=[d22off;d22on];             
             type{i}.fly{j}.rawY2=[zeros(size(d22off,1),1);ones(size(d22on,1),1)];

             type{i}.fly{j}.X=[poff';pon'];
             type{i}.fly{j}.XX=[TPoff';TPon'];
             type{i}.fly{j}.Y=[zeros(size(poff,2),1);ones(size(pon,2),1)];
%             type{i}.fly{j}.X=type{i}.fly{j}.X(:,1:end-1);
%             type{i}.fly{j}.XX=type{i}.fly{j}.XX(:,1:end-1);
                          
             type{i}.fly{j}.X2=[poff2';pon2'];
             type{i}.fly{j}.XX2=[TPoff2';TPon2'];
             type{i}.fly{j}.Y2=[zeros(size(poff2,2),1);ones(size(pon2,2),1)];
%             type{i}.fly{j}.X2=type{i}.fly{j}.X2(:,1:end-1);
%             type{i}.fly{j}.XX2=type{i}.fly{j}.XX2(:,1:end-1);


             rawX{i}=[rawX{i};d2off;d2on];             
             rawY{i}=[rawY{i};zeros(length(idxoff),1);ones(length(idxon),1)];
             
             rawX2{i}=[rawX2{i};d22off;d22on];             
             rawY2{i}=[rawY2{i};zeros(size(d22off,1),1);ones(size(d22on,1),1)];
             
             X{i}=[X{i};poff';pon'];
             XX{i}=[XX{i};TPoff';TPon'];
             Y{i}=[Y{i};zeros(size(poff,2),1);ones(size(pon,2),1)];
             
             X2{i}=[X2{i};poff2';pon2'];
             XX2{i}=[XX2{i};TPoff2';TPon2'];
             Y2{i}=[Y2{i};zeros(size(poff2,2),1);ones(size(pon2,2),1)];
             
             
            [type{i}.fly{j}.rawB,dev,stats] = glmfit(type{i}.fly{j}.rawX,type{i}.fly{j}.rawY,'binomial','link','logit','Options',opts);
            type{i}.fly{j}.rawYhat=glmval(type{i}.fly{j}.rawB,type{i}.fly{j}.rawX,'logit');
%            type{i}.fly{j}.PCraw = sum((type{i}.fly{j}.rawYhat>1/2)==type{i}.fly{j}.rawY)/length(type{i}.fly{j}.rawY);
            type{i}.fly{j}.PCraw = max(sum((type{i}.fly{j}.rawYhat>dc)==repmat(type{i}.fly{j}.rawY,1,length(dc)))/length(type{i}.fly{j}.rawY));

            [type{i}.fly{j}.rawB2,dev,stats] = glmfit(type{i}.fly{j}.rawX2,type{i}.fly{j}.rawY2,'binomial','link','logit','Options',opts);
            type{i}.fly{j}.rawYhat2=glmval(type{i}.fly{j}.rawB2,type{i}.fly{j}.rawX2,'logit');
%            type{i}.fly{j}.PCraw2 = sum((type{i}.fly{j}.rawYhat2>1/2)==type{i}.fly{j}.rawY2)/length(type{i}.fly{j}.rawY2);
            type{i}.fly{j}.PCraw2 = max(sum((type{i}.fly{j}.rawYhat2>dc)==repmat(type{i}.fly{j}.rawY2,1,length(dc)))/length(type{i}.fly{j}.rawY2));

            [type{i}.fly{j}.B,dev,stats] = glmfit(type{i}.fly{j}.X,type{i}.fly{j}.Y,'binomial','link','logit','Options',opts);
            type{i}.fly{j}.Yhat=glmval(type{i}.fly{j}.B,type{i}.fly{j}.X,'logit');
%            type{i}.fly{j}.PC = sum((type{i}.fly{j}.Yhat>1/2)==type{i}.fly{j}.Y)/length(type{i}.fly{j}.Y);
            type{i}.fly{j}.PC = max(sum((type{i}.fly{j}.Yhat>dc)==repmat(type{i}.fly{j}.Y,1,length(dc)))/length(type{i}.fly{j}.Y));

            [type{i}.fly{j}.B2,dev,stats] = glmfit(type{i}.fly{j}.X2,type{i}.fly{j}.Y2,'binomial','link','logit','Options',opts);
            type{i}.fly{j}.Yhat2=glmval(type{i}.fly{j}.B2,type{i}.fly{j}.X2,'logit');
%            type{i}.fly{j}.PC2 = sum((type{i}.fly{j}.Yhat2>1/2)==type{i}.fly{j}.Y2)/length(type{i}.fly{j}.Y2);
            type{i}.fly{j}.PC2 = max(sum((type{i}.fly{j}.Yhat2>dc)==repmat(type{i}.fly{j}.Y2,1,length(dc)))/length(type{i}.fly{j}.Y2));

            [type{i}.fly{j}.TPB,dev,stats] = glmfit(type{i}.fly{j}.XX,type{i}.fly{j}.Y,'binomial','link','logit','Options',opts);
            type{i}.fly{j}.TPYhat=glmval(type{i}.fly{j}.TPB,type{i}.fly{j}.XX,'logit');
%            type{i}.fly{j}.TPPC = sum((type{i}.fly{j}.TPYhat>1/2)==type{i}.fly{j}.Y)/length(type{i}.fly{j}.Y);
            type{i}.fly{j}.TPPC = max(sum((type{i}.fly{j}.TPYhat>dc)==repmat(type{i}.fly{j}.Y,1,length(dc)))/length(type{i}.fly{j}.Y));

            [type{i}.fly{j}.TPB2,dev,stats] = glmfit(type{i}.fly{j}.XX2,type{i}.fly{j}.Y2,'binomial','link','logit','Options',opts);
            type{i}.fly{j}.TPYhat2=glmval(type{i}.fly{j}.TPB2,type{i}.fly{j}.XX2,'logit');
%            type{i}.fly{j}.TPPC2 = sum((type{i}.fly{j}.TPYhat2>1/2)==type{i}.fly{j}.Y)/length(type{i}.fly{j}.Y2);
            type{i}.fly{j}.TPPC2 = max(sum((type{i}.fly{j}.TPYhat2>dc)==repmat(type{i}.fly{j}.Y2,1,length(dc)))/length(type{i}.fly{j}.Y2));

        end
        


        [TPB{i},dev,stats] = glmfit(XX{i},Y{i},'binomial','link','logit','Options',opts);
        TPYhat{i}=glmval(TPB{i},XX{i},'logit');
%        TPpercentCorrect(i)=sum((TPYhat{i}>1/2)==Y{i})/length(Y{i})
        [TPpercentCorrect(i),loc]=max(sum(((TPYhat{i}>=dc) == repmat(Y{i},1,length(dc))))/length(Y{i}));
%        dc(loc)

        [TPB2{i},dev,stats] = glmfit(XX2{i},Y2{i},'binomial','link','logit','Options',opts);
        TPY2hat{i}=glmval(TPB2{i},XX2{i},'logit');
%        B2percentCorrect(i)=sum((TPY2hat{i}>1/2)==Y2{i})/length(Y2{i})
        [TP2percentCorrect(i),loc]=max(sum((TPY2hat{i}>=dc) == repmat(Y2{i},1,length(dc)))/length(Y2{i}));
%        dc(loc)

        [B{i},dev,stats] = glmfit(X{i},Y{i},'binomial','link','logit','Options',opts);
        Yhat{i}=glmval(B{i},X{i},'logit');
%        BpercentCorrect(i)=sum((Yhat{i}>1/2)==Y{i})/length(Y{i})
        [BpercentCorrect(i),loc]=max(sum(((Yhat{i}>=dc) == repmat(Y{i},1,length(dc))))/length(Y{i}));
%        dc(loc)

        [B2{i},dev,stats] = glmfit(X2{i},Y2{i},'binomial','link','logit','Options',opts);
        Y2hat{i}=glmval(B2{i},X2{i},'logit');
%        B2percentCorrect(i)=sum((Y2hat{i}>1/2)==Y2{i})/length(Y2{i})
        [B2percentCorrect(i),loc]=max(sum((Y2hat{i}>=dc) == repmat(Y2{i},1,length(dc)))/length(Y2{i}));
%        dc(loc)
        
        [Braw{i},dev,stats] = glmfit(rawX{i},rawY{i},'binomial','link','logit','Options',opts);
        rawYhat{i}=glmval(Braw{i},rawX{i},'logit');
%        rawpercentCorrect(i)=sum((rawYhat{i}>1/2)==rawY{i})/length(rawY{i})
        [rawpercentCorrect(i),loc]=max(sum(((rawYhat{i}>=dc) == repmat(rawY{i},1,length(dc))))/length(rawY{i}));
%        dc(loc)

        [Braw2{i},dev,stats] = glmfit(rawX2{i},rawY2{i},'binomial','link','logit','Options',opts);
        rawYhat2{i}=glmval(Braw2{i},rawX2{i},'logit');
%        rawpercentCorrect2(i)=sum((rawYhat2{i}>1/2)==rawY2{i})/length(rawY2{i})        
        [rawpercentCorrect2(i),loc]=max(sum(((rawYhat2{i}>=dc) == repmat(rawY2{i},1,length(dc))))/length(rawY2{i}));
%        dc(loc)
        
        for j=1:length(type{i}.fly)
            type{i}.fly{j}.rawYhatALL=glmval(Braw{i},type{i}.fly{j}.rawX,'logit');
            type{i}.fly{j}.PCrawALL = sum((type{i}.fly{j}.rawYhatALL>1/2)==type{i}.fly{j}.rawY)/length(type{i}.fly{j}.rawY);

            type{i}.fly{j}.rawYhat2ALL=glmval(Braw2{i},type{i}.fly{j}.rawX2,'logit');
            type{i}.fly{j}.PCraw2ALL = sum((type{i}.fly{j}.rawYhat2ALL>1/2)==type{i}.fly{j}.rawY2)/length(type{i}.fly{j}.rawY2);

            type{i}.fly{j}.YhatALL=glmval(B{i},type{i}.fly{j}.X,'logit');
            type{i}.fly{j}.PCALL = sum((type{i}.fly{j}.YhatALL>1/2)==type{i}.fly{j}.Y)/length(type{i}.fly{j}.Y);

            type{i}.fly{j}.Yhat2ALL=glmval(B2{i},type{i}.fly{j}.X2,'logit');
            type{i}.fly{j}.PC2ALL = sum((type{i}.fly{j}.Yhat2ALL>1/2)==type{i}.fly{j}.Y2)/length(type{i}.fly{j}.Y2);

            type{i}.fly{j}.TPYhatALL=glmval(TPB{i},type{i}.fly{j}.XX,'logit');
            type{i}.fly{j}.TPPCALL = sum((type{i}.fly{j}.TPYhatALL>1/2)==type{i}.fly{j}.Y)/length(type{i}.fly{j}.Y);

            type{i}.fly{j}.TPYhat2ALL=glmval(TPB2{i},type{i}.fly{j}.XX2,'logit');
            type{i}.fly{j}.TPPC2ALL = sum((type{i}.fly{j}.TPYhat2ALL>1/2)==type{i}.fly{j}.Y2)/length(type{i}.fly{j}.Y2);

            type{i}.fly{j}.chance2 = mean(type{i}.fly{j}.Y2);
            type{i}.fly{j}.chance2 = max(type{i}.fly{j}.chance2,1-type{i}.fly{j}.chance2);            
            type{i}.fly{j}.chance = mean(type{i}.fly{j}.Y);
            type{i}.fly{j}.chance = max(type{i}.fly{j}.chance,1-type{i}.fly{j}.chance);
        end
    end
end
chance = 1-[mean(Y{1}),mean(Y{2}),mean(Y{3}),mean(Y{4})]

figure
kk=gcf;
kk=kk.Number-1;
cc=prism(4);
for n=1:length(type)
for i=1:length(type{n}.fly)
    PCraw{n}(i)=type{n}.fly{i}.PCraw/type{n}.fly{i}.chance;
    PC{n}(i)=type{n}.fly{i}.PC/type{n}.fly{i}.chance;
    PCraw2{n}(i)=type{n}.fly{i}.PCraw2/type{n}.fly{i}.chance2;
    PC2{n}(i)=type{n}.fly{i}.PC2/type{n}.fly{i}.chance2;
    
    TPPC{n}(i)=type{n}.fly{i}.TPPC/type{n}.fly{i}.chance;
    TPPC2{n}(i)=type{n}.fly{i}.TPPC2/type{n}.fly{i}.chance2;
    TPPCALL{n}(i)=type{n}.fly{i}.TPPCALL/type{n}.fly{i}.chance;
    TPPC2ALL{n}(i)=type{n}.fly{i}.TPPC2ALL/type{n}.fly{i}.chance2;
    
end
figure(kk+1);scatter(PCraw{n},PC{n},30*ones(size(PC{n})),repmat(cc(n,:),length(PC{n}),1));hold on
figure(kk+2);scatter(PCraw2{n},PC2{n},30*ones(size(PC2{n})),repmat(cc(n,:),length(PC2{n}),1));hold on
figure(kk+3);scatter(PCraw{n},TPPC{n},30*ones(size(PC{n})),repmat(cc(n,:),length(PC{n}),1));hold on
figure(kk+4);scatter(PCraw2{n},TPPC2{n},30*ones(size(PC2{n})),repmat(cc(n,:),length(PC2{n}),1));hold on

end

figure(kk+1);
refline(1)
title('Single Fly Fits')
xlabel('PCraw')
ylabel('PC')
legend('1','2','3','4')
hold off
figure(kk+2)
refline(1)
title('Single Fly Fits')
xlabel('PCraw2')
ylabel('PC2')
legend('1','2','3','4')
hold off;
figure(kk+3);
refline(1)
title('Single Fly Fits')
xlabel('PCraw')
ylabel('TPPC')
legend('1','2','3','4')
hold off
figure(kk+4)
refline(1)
title('Single Fly Fits')
xlabel('PCraw2')
ylabel('TCPC2')
legend('1','2','3','4')
hold off;


for n=1:length(type)
for i=1:length(type{n}.fly)
    PCrawALL{n}(i)=type{n}.fly{i}.PCrawALL/type{n}.fly{i}.chance;
    PCALL{n}(i)=type{n}.fly{i}.PCALL/type{n}.fly{i}.chance;
    PCraw2ALL{n}(i)=type{n}.fly{i}.PCraw2ALL/type{n}.fly{i}.chance;
    PC2ALL{n}(i)=type{n}.fly{i}.PC2ALL/type{n}.fly{i}.chance;
end
figure(kk+5);scatter(PCrawALL{n},PCALL{n},30*ones(size(PCALL{n})),repmat(cc(n,:),length(PCALL{n}),1));hold on
figure(kk+6);scatter(PCraw2ALL{n},PC2ALL{n},30*ones(size(PC2ALL{n})),repmat(cc(n,:),length(PC2ALL{n}),1));hold on
figure(kk+7);scatter(PCrawALL{n},TPPCALL{n},30*ones(size(PCALL{n})),repmat(cc(n,:),length(PCALL{n}),1));hold on
figure(kk+8);scatter(PCraw2ALL{n},TPPC2ALL{n},30*ones(size(PC2ALL{n})),repmat(cc(n,:),length(PC2ALL{n}),1));hold on

end

figure(kk+5);
refline(1)
title('Cluster Fly Fits')
xlabel('PCraw')
ylabel('PC')
legend('1','2','3','4')
hold off
figure(kk+6)
refline(1)
title('Cluster Fly Fits')
xlabel('PCraw2')
ylabel('PC2')
legend('1','2','3','4')
hold off;

figure(kk+7);
refline(1)
title('Cluster Fly Fits')
xlabel('PCraw')
ylabel('TPPC')
legend('1','2','3','4')
hold off
figure(kk+8)
refline(1)
title('Cluster Fly Fits')
xlabel('PCraw2')
ylabel('TPPC2')
legend('1','2','3','4')
hold off;






% 
% 
% [m,idx]=max(model.p);
% cc=jet(model.Qdim);
% plotidx=0;
% for i=1:model.NC
%     idx2=find(idx==i);  % idx of datapoints assigned to cluster i
%     np=length(idx2);
%     if(np>0)
%         plotidx=plotidx+1;
%         
%         px=ceil(sqrt(np));
%         py=ceil(np/px);
%         for j=1:np
%             d1=data{idx2(j)};
% %            o1=odoron{idx2(j)};
%             o1=ones(size(idx2(j)));
%             
%             p1=model.HHMMs{i}.p{idx2(j)};
%             clear p2
%             for k=1:model.Qdim
%                 p2(k,:)=sum(p1(model.HHMMs{i}.Aidx{k},:),1);
%             end
%             [m2,idx3]=max(p2);
%             figure(2*plotidx-1);        
%             subplot(px,py,j), scatter(d1(1,:),d1(2,:),3*ones(size(d1(1,:))),cc(idx3,:))
%             subplot(px,py,j), title(['Cluster ',num2str(i), 'with P = ',num2str(model.p(i,idx2(j)))])
%              figure(2*plotidx)
%              idxoff=find(o1==0);
%              idxon=find(o1==1);
%              idxon=idxon(1:min(1000,length(idxon)));
%              subplot(np,3,3*j-2), bar(mean(p2(:,idxoff)')), title('odor off')
%              subplot(np,3,3*j-1), bar(mean(p2(:,idxon)')), title('odor on')
%              subplot(np,3,3*j),   bar(mean(p2(:,idxon)')-mean(p2(:,idxoff)')), title('difference')
%              
%         end
%     end
% end
% 
% 
% 
