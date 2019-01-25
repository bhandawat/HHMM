
load('C:\Users\brain_000\Dropbox\Penalty Kick Ass\data - 8.8.15.mat')
%first 2000 are musimol
% data = data(randi(length(data),1,400));
NC=50;
dim=6;

D=6; obsTypes{1}.idx=[1,2,3,5,6,7];
%D=6; obsTypes{1}.idx=[1,2,3];
obsTypes{1}.dist='mvn';

model = mixHMM(NC,dim,D,obsTypes,ones(NC,1)/NC*4,ones(dim,dim)/dim,ones(dim,1)/dim);
%model1 = HMM(dim*NC,D,obsTypes,ones(dim*NC,dim*NC),ones(dim*NC,1));
maxiters = 200;
%could also use kmeans applied to last N data points to preprocess to get
%good initial assigments.

k=0;
tic
hours=18;
model.initialize(data,40);
%model1.initialize(data);

tic
%model1.update(data)
model.update(data,1)
['1 iteration completed in ',num2str(toc), 'seconds']
model.prune
while(k<maxiters & toc<hours*60*60)
    k=k+1  
 %   model1.update(data)
    model.update(data,1)
    model.pi.mean'
    if(mod(k,10)==0)
        model.prune;
    end
    model.plotclusters(data,1)
    
    figure(2)
    subplot(2,1,1), bar(sum(model.p(:,1:2000,:),2))
    title(['ELBO = ',num2str(model.L),' at iter = ',num2str(k)])
    subplot(2,1,2), bar(sum(model.p(:,2001:end),2))
    subplot(2,1,1), title(['iter = ',num2str(k)])

    hmusc = sum(model.p(:,1:2000),2)/2000;
    idx=find(hmusc>0);
    ENTmusc= -sum(hmusc(idx).*log2(hmusc(idx)))

    hsaline = sum(model.p(:,2001:end),2)/2000;
    idx = find(hsaline>0);
    ENTsaline= -sum(hsaline(idx).*log2(hsaline(idx)))
    
end


