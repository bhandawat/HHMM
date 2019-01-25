for i=1:34
    a(i)=max(data{1,i}(16,:));
end
b=max(a);
factor=20;
for i=1:34
  data{1,i}(16,:)=data{1,i}(16,:)/(b*1/factor);  
end
filename=['WT_ACV0_May2017_' num2str(factor) '.mat'];
save(filename, 'data');