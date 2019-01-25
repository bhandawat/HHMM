function [] = plotColorSteps(data,x,cc)

if isempty(cc)
    cc=varycolor(max(data)-min(data)+1);
end

if isempty(x)
    x = 1:1:length(data);
end
assert(numel(x)==numel(data),'x and y have inconsistent number of values');

ts = [1, find(abs(diff(data))>0)+1, length(data)];
if(ts(end)==ts(end-1))
    ts(end)=[];
end
if(ts(1)==ts(2))
    ts(1)=[];
end

hold on
ylim([min(data)-1 max(data)+1])
for k=1:length(ts)-1
    if data(ts(k))==0                % plotting vertical lines between state transitions white
        plot(x(ts(k):ts(k+1)-1),data(ts(k):ts(k+1)-1),'Color','w');
    else                                                            % plotting instances where state transitions don't occur (between ts(t) and ts(t+1))
        plot(x(ts(k):ts(k+1)-1),data(ts(k):ts(k+1)-1),'Color',cc(data(ts(k)),:),'Linewidth',3);
        plot(x(ts(k+1)-1:ts(k+1)),data(ts(k):ts(k)+1),'Color',cc(data(ts(k)),:),'Linewidth',3);
    end
end
end