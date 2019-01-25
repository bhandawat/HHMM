function [] = HLSObsModContour2(data,model,likely_state_by_fly,HMMsortNdx,clusterNdx)
%%
conf = 0.85; % confidence bound for computing contour for empirical data from the model.
[~,cidx]=max(model.p);
obsdims=model.obsTypes{1}.idx;
ct2 = 0;

for modNA=1:1
    
    if model.NA(modNA)>3
        i = modNA;
        idx=find(cidx==i);
        if(isempty(idx))
        else
            m1=model.HMMs{i};
            Qbar=zeros(m1.dim,1);
            for k=1:m1.dim
                emp{i,k}=[];
            end
            for j=idx
                dtemp=data{j}(obsdims,:);
                stemp = likely_state_by_fly(j,:);
                Qbar=Qbar+sum(m1.p{j},2);                                          % % chance for each HLS over all tracks 10x1 mat
                %sltemp=cell(1,m1.dim);
                for k=1:m1.dim
                    emp{i,k}=[emp{i,k},dtemp(:,(stemp==k))];
                    %sltempAll = find(stemp==k);
                    
                end
            end
            Qbar=Qbar/sum(Qbar);                                                    % distribution of HLS for entire track
            
            figure(i);set(gcf,'position',[9 49 824 918])
            kk=0;
            pHL = cell(8,1);
            for k=HMMsortNdx
                if(Qbar(k)<0.0001)
                    
                else
                    
                    xmin=Inf;xmax=-Inf;
                    ymin=Inf;ymax=-Inf;
                    
                    mu{k}=m1.obsModels{1,k}.dists{1}.mu;
                    Sigma{k}=m1.obsModels{1,k}.dists{1}.E.Sigma;
                    invSigma{k}=inv(Sigma{k});
                    
                    xmin=min(xmin,mu{k}(1)-4*sqrt(Sigma{k}(1,1)));
                    xmax=max(xmax,mu{k}(1)+4*sqrt(Sigma{k}(1,1)));
                    ymin=min(ymin,mu{k}(2)-4*sqrt(Sigma{k}(2,2)));
                    ymax=max(ymax,mu{k}(2)+4*sqrt(Sigma{k}(2,2)));
                    
                    X{k}=round(xmin,2)-0.1:0.01:round(xmax,2)+0.1;
                    Y{k}=round(ymin,2)-0.1:0.01:round(ymax,2)+0.1;
                    
                    x=X{k}'*ones(1,length(Y{k}));
                    y=ones(length(X{k}),1)*Y{k};
                    p=zeros(size(x));
                    
                    lnp = -1/2*( (x-mu{k}(1)).^2*invSigma{k}(1,1) + (y-mu{k}(2)).^2*invSigma{k}(2,2) + 2*(x-mu{k}(1)).*(y-mu{k}(2))*invSigma{k}(1,2));
                    lnp = lnp-1/2*log(det(2*pi*Sigma{k}));
                    p=p+Qbar(k)*exp(lnp);
                    prepad = round([(8-X{k}(end))./0.01, (4.5-Y{k}(end))./0.01]);
                    postpad = round([(X{k}(1)+4)./0.01, (Y{k}(1)+4.5)./0.01]);
                    pLL{k} = padarray(p,postpad,0,'pre');
                    pLL{k} = padarray(pLL{k},prepad,0,'post');
                    
                    
                end
            end
            pLLSorted = pLL(HMMsortNdx);
            
            figure;
            x=-4:0.01:8;y = -4.5:0.01:4.5;
            c=varycolor(10);
            for HLS = 1:10
                pHLS = zeros(size(pLLSorted{1}));
                for LLS = 1:length(clusterNdx{HLS})
                    pHLS = pHLS+pLLSorted{clusterNdx{HLS}(LLS)};
                end

                pHLS=pHLS/sum(pHLS(:));                                                  % normalize
                m=sort(pHLS(:),'descend');
                v=max(find(cumsum(m)<conf));
                v=m(v);
                contour(x,y,pHLS',[v v],'LineWidth',2,'Color',c(HLS,:));hold on
                axis([-1,5,-1.5,1.5])

            end
            title(['HL States (' num2str(HLS) ')'])
            set(gcf,'position',[9 49 824 918])
        end
    end
end
end