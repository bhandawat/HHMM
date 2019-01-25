function [RiV,RV,I,C, vMinWeights] = iVat(R,VATflag,cutflag)
% Example function call: [RiV] = iVat(RV);
%
% *** Input Parameters ***
% @param R (n*n double): dissimilarity data input
% @param R (n*D double): vector input (R is converted to sq. Euclidean
% distance)
% @param VATflag (boolean): TRUE - R is VAT-reordered
% @param cutflag (boolean): TRUE - return the sinlge-linkage partition
%                           according to chosen cut (use mouse to choose
%                           the cut distance value from the visualization)
%
% *** Output Values ***
% @value RV (n*n double): VAT-reordered dissimilarity data
% @value RiV (n*n double): iVAT-transformed dissimilarity data
% @value I (n int): reordering indices
% @value C (n int): cluster indicator vector (cutflag == 1); note: C is a
%                   single-linkage partition
%
% Modified by Jeffrey Chan, 6/2013
%

if(nargin==1)
    VATflag = 0;
    cutflag = 0;
elseif(nargin==2)
    cutflag = 0;
end;
if(cutflag)
    VATflag=1;
end;

[N,M]=size(R);
if(N~=M)
    R=distance2(R,R);
end;

if(VATflag),
    RV=R;
    RiV=zeros(N);
    for r=2:N,
        c=1:r-1;
        [y i]=min(RV(r,1:r-1));
        RiV(r,c)=y;
        cnei=c(c~=i);
        RiV(r,cnei)=max([RiV(r,cnei); RiV(i,cnei)]);
        vMinWeights = RiV(r,cnei);
        RiV(c,r)=RiV(r,c)';
    end;
else
%     [RV,C,I,~,~]=VAT(R);
    [RV,I,C]=Vat(R);
    
    RiV=zeros(N);
    for r=2:N,
        c=1:r-1;
        RiV(r,c)=RV(r,C(r));
        cnei=c(c~=C(r));
        RiV(r,cnei)=max([RiV(r,cnei); RiV(C(r),cnei)]);
        vMinWeights = RiV(r,cnei);
        RiV(c,r)=RiV(r,c)';
    end;
end;

%        figure;
%        colormap(gray);
%        imagesc(RiV);


if(cutflag)
    notdone=1; cmin = min(min(triu(RiV,1))); cmax=max(max(triu(RiV,1)));
    while(notdone)
        figure; set(gcf,'DefaultLineLineWidth',2,'DefaultLineColor','y');
        imagesc(RiV); axis image; axis([-0.5 N+1 -0.5 N+1]); axis off; colormap gray; caxis([cmin cmax]);
        disp('Use mouse to click on the cut value (the darkest value between blocks)');
        [x,y]=ginput(1);
        cutval = RiV(round(y),round(x));
        i = find(diag(RiV,1)>=cutval)+0.5;
        line('XData',[0.5 0.5],'YData',[0.5 i(1)],'LineWidth',4,'Color','g');
        line('XData',[0.5 i(1)],'YData',[0.5 0.5],'LineWidth',4,'Color','g');
        line('XData',[i(1) i(1)],'YData',[0.5 i(1)],'LineWidth',4,'Color','g');
        line('XData',[0.5 i(1)],'YData',[i(1) i(1)],'LineWidth',4,'Color','g');
        for j=1:length(i)-1,
            line('XData',[i(j) i(j)],'YData',[i(j) i(j+1)],'LineWidth',4,'Color','g');
            line('XData',[i(j) i(j+1)],'YData',[i(j) i(j)],'LineWidth',4,'Color','g');
            line('XData',[i(j+1) i(j+1)],'YData',[i(j) i(j+1)],'LineWidth',4,'Color','g');
            line('XData',[i(j) i(j+1)],'YData',[i(j+1) i(j+1)],'LineWidth',4,'Color','g');
        end;
        line('XData',[i(end) i(end)],'YData',[i(end) N+0.5],'LineWidth',4,'Color','g');
        line('XData',[i(end) N+0.5],'YData',[i(end) i(end)],'LineWidth',4,'Color','g');
        line('XData',[N+0.5 N+0.5],'YData',[i(end) N+0.5],'LineWidth',4,'Color','g');
        line('XData',[i(end) N+0.5],'YData',[N+0.5 N+0.5],'LineWidth',4,'Color','g');
        
        notdone=input('Would you like the select again? (0 - no, 1 - yes) ');
    end;
    ind = 1:N; C=zeros(N,1);
    C(ind<i(1))=1;
    for j=2:length(i),
        C(ind<i(j) & ind>i(j-1))=j;
    end;
    C(ind>i(end))=j+1;
else
    C=[];
end;