function [] = convertToExcel()

for i = 1:34
    x = data{i}(1,:)';
    y = data{i}(2,:)';
    speed = data{i}(7,:)';
    curvature = data{i}(8,:)';
    vParallel = data{i}(15,:)'./2;
    vPerpendicular = data{i}(16,:)';

    %T = table(x,y,speed,curvature,vParallel,vPerpendicular);
    T = table(x,y);
    filename = ['fly' num2str(i) '_Data.xlsx'];
    writetable(T,filename,'Sheet',1,'Range','A1')
end

speed = [];speed2 = [];speed3 = [];
for i = 1:34
    x = data{i}(1,:)';
    y = data{i}(2,:)';
    speed(:,i) = sqrt(diff(x.*3.2).^2+diff(y.*3.2).^2)*30;
    speed2(:,i) = sqrt((data{i}(15,:)'/2).^2+data{i}(16,:)'.^2);
    speed3(:,i) = data{i}(7,:)';
end

figure;
subplot(3,1,1);histogram(speed(:),[0:0.1:3]);
subplot(3,1,2);histogram(speed2(:),[0:0.1:3]);
subplot(3,1,3);histogram(speed3(:)./1.8750,[0:0.1:3]);



end
