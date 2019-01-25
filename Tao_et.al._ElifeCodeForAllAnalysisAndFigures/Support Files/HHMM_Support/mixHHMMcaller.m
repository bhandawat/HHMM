NC=length(data);
Qdim=8;
dim=4;
D=2;
obsTypes{1}.idx=[15,16];
obsTypes{1}.dist='mvn';
%obsTypes{2}.idx=9;
%obsTypes{2}.dist='gamma';
% obsTypes{1}.idx=13;
% obsTypes{1}.dist='gamma';
% obsTypes{2}.idx=14;
% obsTypes{2}.dist='normal';

model = mixHHMM(NC,Qdim,dim,D,obsTypes,ones(NC,1),eye(Qdim)+ones(Qdim)/Qdim, ones(Qdim,1)/Qdim,eye(dim)+ones(dim)/dim,ones(dim,1)/dim);
maxiters = 20;
%could also use kmeans applied to last N data points to preprocess to get
%good initial assigments.

k=0;
hours=15;
DL=Inf;

model.initialize(data,200)
model.update(data,5,4)
model.prune();
model.update(data,1)
model.plotclusters(data,1);

%while(k<maxiters & toc<hours*60*60 & DL > 1)
maxiters=20; 
tic
while(k<maxiters)
    k=k+1
    if(k<=maxiters/2)
        model.prune();
        toc
        model.update(data,10,4);
        toc
        model.split(data,50);
        toc
        model.update(data,10,4);
        toc
        model.fillunused(data,50);
        toc
        model.update(data,10,4);
        toc
        DL=model.update(data,1)
        toc
    else
        DL=model.update(data,10,4)
    end
        model.plotclusters(data,1)
end
%%
[m,idx]=max(model.p);
cc=jet(model.Qdim);
plotidx=0;
for i=1:model.NC
    idx2=find(idx==i);  % idx of datapoints assigned to cluster i
    np=length(idx2);
    if(np>0)
        plotidx=plotidx+1;
        
        px=ceil(sqrt(np));
        py=ceil(np/px);
        for j=1:np
            d1=data{idx2(j)};
%            o1=odoron{idx2(j)};
            o1=ones(size(idx2(j)));
            
            p1=model.HHMMs{i}.p{idx2(j)};
            clear p2
            for k=1:model.Qdim
                p2(k,:)=sum(p1(model.HHMMs{i}.Aidx{k},:),1);
            end
            [m2,idx3]=max(p2);
            figure(2*plotidx-1);        
            subplot(px,py,j), scatter(d1(1,:),d1(2,:),3*ones(size(d1(1,:))),cc(idx3,:))
            subplot(px,py,j), title(['Cluster ',num2str(i), 'with P = ',num2str(model.p(i,idx2(j)))])
             figure(2*plotidx)
             idxoff=find(o1==0);
             idxon=find(o1==1);
             idxon=idxon(1:min(1000,length(idxon)));
             subplot(np,3,3*j-2), bar(mean(p2(:,idxoff)')), title('odor off')
             subplot(np,3,3*j-1), bar(mean(p2(:,idxon)')), title('odor on')
             subplot(np,3,3*j),   bar(mean(p2(:,idxon)')-mean(p2(:,idxoff)')), title('difference')
             
        end
    end
end



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
