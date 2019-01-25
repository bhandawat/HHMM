function DiffusionAnalysis(printFig)
%DiffusionAnalysis Diffusion Analysis is a rudimentary method to compare
%the quality of models. In this code, empirical tracks are generated using
%fly trajectory data. Synthetic tracks are generated using a Run and Tumble
%model and a fitted HHMM model.

% Inputs:
%    printFig: 1 for print to ps, otherwise no

% Outputs:
%    N/A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure List:
%       1.) Sample autocorrelation for speed and curvature
%       2.) Sample and net HLS probability from generative HHMM and fitted
%       empirical HHMM tracks.
%       3-4.) 12 sample empirical and synthetic tracks
%       5.) Diffusion analysis at intervals of 1,2,5,10,20,50,100,300,900
%       frames
%       6.) Jenson-Shannon and Euclidean distances between diffusion
%       analysis distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: 1.) color scheme: red = emp, green = Run and Tumble, blue = HHMM 
%       2.) change scenario to switch between before and during
%       2.) have to change figure labels manually for before and during -
%       change later
%
% 2018, Liangyu Tao

close all
scenario = 1;% 1 = before, other = during

% obtain empirical data
[dataEmp,model,numFliesEmp,first_entry] = loadData('emp');

for i = 1:34
    first_entry(i) = find(dataEmp{i}(11,:)&dataEmp{i}(12,:),1);
end


% generate synthetic and empirical tracks
%[Emp,EmpSC] = EmpAnalysis(dataEmp,numFliesEmp,first_entry,scenario);
[HHMM] = HHMMAnalysis(dataEmp,model,numFliesEmp,first_entry,scenario,printFig);

c = varycolor(4);
% plot model comparisons
figure(3)
[~,~,dJS,dEu] = plotModelComparisons(EmpSC.dist,RT.dist,HHMM.reg.dist,HHMM.reg.dist,c);
% subplot(3,3,1)
% legend('EmpSC','RT','HHMM','HHMM')
% set(gcf,'Position',[849 49 824 918])
% suptitle('Before,0.5')
% if printFig == 1
%     print('-dpsc2',[['PDF_Files/ModelComparison'] '.ps'],'-loose','-append');
% end

% plot Jenson-Shannon Divergence
dist2Consider = [1,2,5,10,20,90,100,300,900];
%dist2Consider = [1,2,3,15:15:900];
figure(4)
subplot(2,1,1);hold on;
for i = 1:2
    plot(dist2Consider,dJS(i,:),'Color',c(i,:))
end
xlabel('Time Offset');ylabel('Jenson-Shannon Divergence')
legend('EmpSC-RT','EmpSC-HHMM')

% plot Euclidean Distance
subplot(2,1,2);hold on;
for i = 1:2
    plot(dist2Consider,dEu(i,:),'Color',c(i,:))
end
xlabel('Time Offset');ylabel('Euclidean Distance (50D)')
legend('EmpSC-RT','EmpSC-HHMM')
set(gcf,'Position',[849 49 824 918])
suptitle('Before,0.5')
if printFig == 1
    print('-dpsc2',[['PDF_Files/ModelComparison'] '.ps'],'-loose','-append');
end
print('-dpdf',[['PDF_Files/ModelComparisonEvery15'] '.pdf']);

plotSampTracks(EmpSC,RT,HHMM.reg,HHMM.reg,c,printFig);

end

function [data,model,numFlies,first_entry] = loadData(type)

if strcmpi(type,'emp')
    % find and load data files
    fnames = dir('Data/*.mat');
    filename = fnames(1).name;
    folder_data = 'Data/';
    load([folder_data filename])
    [data, first_entry] = post_hoc_preprocess_LT(data);
    numFlies = length(data);
end

end

function [Emp,EmpSC] = EmpAnalysis(data,numFlies,first_entry,scenario)

x = zeros(numFlies,size(data{1},2));y = zeros(numFlies,size(data{1},2));
speed = zeros(numFlies,size(data{1},2));curv = zeros(numFlies,size(data{1},2));
for i = 1:numFlies
    x(i,:) = data{i}(3,:);
    y(i,:) = data{i}(4,:);
    speed(i,:) = data{i}(7,:);
    curv(i,:) = data{i}(8,:);
end

xPosAll = cell(34,1);yPosAll = cell(34,1);spdMat = cell(34,1);curvMat = cell(34,1);
if scenario == 1
    % Before case
    distAway = cell(34,1);
    for i = 1:numFlies
        kk = 1;
        for j = 1:300:first_entry(i)-1000
            xTemp = x(i,j:j+1000);
            yTemp = y(i,j:j+1000);
            spdMat{i}(kk,:) = speed(i,j:j+1000);
            curvMat{i}(kk,:) = curv(i,j:j+1000);
            distAway{i}(kk,:) = sqrt((xTemp-xTemp(1)).^2+(yTemp-yTemp(1)).^2);
            xPosAll{i}(kk,:) = xTemp;
            yPosAll{i}(kk,:) = yTemp;
            kk = kk+1;
        end
    end
    
else
    % During case
    distAway = cell(34,1);
    for i = 1:numFlies
        kk = 1;
        for j = first_entry(i):300:length(x)-1001
            xTemp = x(i,j:j+1000);
            yTemp = y(i,j:j+1000);
            spdMat{i}(kk,:) = speed(i,j:j+1000);
            curvMat{i}(kk,:) = curv(i,j:j+1000);
            distAway{i}(kk,:) = sqrt((xTemp-xTemp(1)).^2+(yTemp-yTemp(1)).^2);
            xPosAll{i}(kk,:) = xTemp;
            yPosAll{i}(kk,:) = yTemp;
            kk = kk+1;
        end
    end
end

spdMat = cell2mat(spdMat);
curvMat = cell2mat(curvMat);
xPosAll = cell2mat(xPosAll);
yPosAll = cell2mat(yPosAll);
distAll = cell2mat(distAway);

[xPosAllSyn,yPosAllSyn,distAllSyn] = XYPosCalc(xPosAll(:,1),yPosAll(:,1),spdMat,curvMat);
EmpSC.x = xPosAllSyn;EmpSC.y = yPosAllSyn;EmpSC.dist = distAllSyn;
Emp.x = xPosAll;Emp.y = yPosAll;Emp.dist = distAll;

%plot autocorrelation
figure(100)
acf = autocorr(spdMat(1,:),500);
subplot(2,1,1);hold on;plot(acf,'r')
xlabel('frames');ylabel('Sample autocorrelation');title('Speed')
acf = autocorr(curvMat(1,:),500);
subplot(2,1,2);hold on;plot(acf,'r')
xlabel('frames');ylabel('Sample autocorrelation');title('Curvature')


end

% Run and Tumble Model
function [RT] = RTanalysis(data,numFlies,first_entry,scenario)

x = zeros(numFlies,size(data{1},2));
for i = 1:numFlies
    x(i,:) = data{i}(3,:);
    y(i,:) = data{i}(4,:);
end

if scenario == 1
    % Before case
    xTemp = [];yTemp = [];
    for i = 1:numFlies
        j = 1:300:first_entry(i)-1000;
        xTemp = [xTemp; x(i,j)'];
        yTemp = [yTemp; y(i,j)'];
    end
else
    % Before case
    xTemp = [];yTemp = [];
    for i = 1:numFlies
        j = first_entry(i):300:length(x)-1001;
        xTemp = [xTemp; x(i,j)'];
        yTemp = [yTemp; y(i,j)'];
    end
end

cutoffCurv = 120./30.*pi./180;
flies = 1:34;

[curvDistBefore,curvDistDuring,distDistBefore,distDistDuring,durDistBefore,durDistDuring] = durCurvProcess(data,flies,cutoffCurv);

if scenario == 1
    % Before case
    cutOff = 0.5;
    [synthFlyAll] = RunAndTumble(distDistBefore,curvDistBefore,durDistBefore,length(xTemp),cutOff,1,1000,xTemp.*3.2,yTemp.*3.2);
    distAll = sqrt((synthFlyAll.x-synthFlyAll.x(:,1)).^2+(synthFlyAll.y-synthFlyAll.y(:,1)).^2)./3.2;
else
    % During case
    cutOff = 0.75;
    [synthFlyAll] = RunAndTumble(distDistDuring,curvDistDuring,durDistDuring,length(xTemp),cutOff,1,1000,xTemp.*3.2,yTemp.*3.2);
    distAll = sqrt((synthFlyAll.x-synthFlyAll.x(:,1)).^2+(synthFlyAll.y-synthFlyAll.y(:,1)).^2)./3.2;
end
xPosAll = synthFlyAll.x./3.2;
yPosAll = synthFlyAll.y./3.2;

RT.x = xPosAll;RT.y = yPosAll;RT.dist = distAll;

end

function [curvDistBefore,curvDistDuring,distDistBefore,distDistDuring,durDistBefore,durDistDuring] = durCurvProcess(data,flies,cutoffCurv)

%obtain all possible speed and curvatures (distribution)
speed = [];curv = [];
for i = flies
    speed = [speed, data{i}(7,:)];
    curv = [curv, data{i}(8,:)];
end

%obtain all speed and curv before first entrance
speedBefore = [];curvBefore = []; speedAfter = []; curvAfter = [];R_b = [];R_a = [];
for i = flies
    b=find(data{i}(11,1:end-1)==1 & data{i}(12,1:end-1)==1);       % odor on and fly inside
    speedBefore = [speedBefore, data{i}(7,1:b(1)-1)];
    curvBefore = [curvBefore, data{i}(8,1:b(1)-1)];
    speedAfter = [speedAfter, data{i}(7,b(1):end)];
    curvAfter = [curvAfter, data{i}(8,b(1):end)];
    R_b = [R_b, sqrt(data{i}(3,1:b(1)-1).^2+data{i}(4,1:b(1)-1).^2)];
    R_a = [R_a, sqrt(data{i}(3,b(1):end).^2+data{i}(4,b(1):end).^2)];
end

%Before case**************************************************************
[distDistBefore,durDistBefore,curvDistBefore] = LenDurProcessing(speedBefore,curvBefore,cutoffCurv,R_b);
%During case**************************************************************
[distDistDuring,durDistDuring,curvDistDuring] = LenDurProcessing(speedAfter,curvAfter,cutoffCurv,R_a);

end

function [distDist,durDist,curvDist] = LenDurProcessing(speed,curv,cutoffCurv,R)

%finding turns
curv(speed<0.1) = [];           % remove stops designated as <1mm/s
speed(speed<0.1) = [];
[~,locs] = findpeaks(abs(curv),'MinPeakHeight',cutoffCurv);

%removing turns where abs diff in angle is <20 deg
closeAngs = locs'+repmat([-2:1:2],length(locs),1);
diffAng = sum(curv(closeAngs),2);
locs2 = locs(abs(diffAng)>=20.*pi./180);

%removing turns where fly is stationary (speed<1mm/s)
spdAtTurn = speed(locs2);
locs2(10.*spdAtTurn<1) = [];

% %removing turns where fly is at border (r>2.5 mm from border edge (1))
% distAtTurn = R(locs2);
% locs2(distAtTurn.*3.2>3.2-0.25) = [];
 curvDist = curv(locs2);

%calculate length distribution
distDist = zeros(1,length(locs2)-1);
R_d = abs(diff(R));
for i = 1:length(locs2)-1
    distDist(i) = sum(abs(speed(locs2(i):locs2(i+1)-1))./30);
    %distDist(i) = sum(R_d(locs2(i):locs2(i+1)-1)).*3.2;
end
durDist = diff(locs2);
curvDist = curvDist(1:end-1);

end

% HHMM Model
function [HHMM] = HHMMAnalysis(data,model,numFlies,first_entry,scenario,printFig)
load('F:\Liangyu Tao\HHMM\Final Code\mat_Files\HLS_by_fly.mat')
load('TP.mat')
load('LLSProb.mat')
load('LLSSCemp.mat')

x = zeros(numFlies,size(data{1},2));
for i = 1:numFlies
    x(i,:) = data{i}(3,:);
    y(i,:) = data{i}(4,:);
end

if scenario == 1
    % Before case
    xTemp = [];yTemp = [];initState = [];allState = [];
    for i = 1:numFlies
        j = 1:300:first_entry(i)-1000;
        xTemp = [xTemp; x(i,j)'];
        yTemp = [yTemp; y(i,j)'];
        initState = [initState; likely_high_state_sorted(i,j)'];
        for k = 1:length(j)
            allState = [allState; likely_high_state_sorted(i,j(k):j(k)+1000)];
        end
    end
else
    % During case
    xTemp = [];yTemp = [];initState = [];allState = [];
    for i = 1:numFlies
        j = first_entry(i):300:length(x)-1001;
        xTemp = [xTemp; x(i,j)'];
        yTemp = [yTemp; y(i,j)'];
        initState = [initState; likely_high_state_sorted(i,j)'];
        for k = 1:length(j)
            allState = [allState; likely_high_state_sorted(i,j(k):j(k)+1000)];
        end
    end
end
initState(initState == 0) = 1;
allState(allState == 0) = 1;

TPAll = zeros(model.Qdim);
for i = 1:34
    TP = model.HHMMs{1}.getQxi(data(i));
    TP = TP{1};
    if scenario == 1
        TPAll = TPAll+sum(TP(:,:,1:first_entry(i)),3);
    else
        TPAll = TPAll+sum(TP(:,:,first_entry(i):end),3);
    end
end
ndx = [9,10,5,3,1,4,2,6,7,8];% hard coded from preprocess
TP_sorted=TPAll(ndx,ndx);
TPAll = TP_sorted./repmat(sum(TP_sorted,1),model.Qdim,1);

[state,stateLL] = syntheticFlyGeneration(initState,TPAll,TP_LL_normalized,LLSProbDist);
[conSpdSpdMat,conSpdCurvMat] = synSCGenRTStyle(state,stateLL,empLLSP_sorted,'constSpeed');
[conCurvSpdMat,conCurvCurvMat] = synSCGenRTStyle(state,stateLL,empLLSP_sorted,'constCurv');
[spdMat,curvMat] = synSpdCurvGen(state,stateLL,empLLSP_sorted);

% plot autocorrelation and HLS distribution
HHMMPlottingFunction(spdMat,curvMat,allState,state,printFig);

[xPosConSpd,yPosConSpd,distConSpd] = XYPosCalc(xTemp,yTemp,conSpdSpdMat,conSpdCurvMat);
[xPosConCurv,yPosConCurv,distConCurv] = XYPosCalc(xTemp,yTemp,conCurvSpdMat,conCurvCurvMat);
[xPosAll,yPosAll,distAll] = XYPosCalc(xTemp,yTemp,spdMat,curvMat);
HHMM = [];
HHMM.conSpd.x = xPosConSpd;HHMM.conSpd.y = yPosConSpd;HHMM.conSpd.dist = distConSpd;
HHMM.conCurv.x = xPosConCurv;HHMM.conCurv.y = yPosConCurv;HHMM.conCurv.dist = distConCurv;
HHMM.reg.x = xPosAll;HHMM.reg.y = yPosAll;HHMM.reg.dist = distAll;

end

function [state,stateLL] = syntheticFlyGeneration(initState,TP_normalized,TP_LL_normalized,LLSProbDist)
%generate HLS distributions
HLSDist = cell(1,10);
for i = 1:10
    for j = 1:10
        HLSDist{i} = [HLSDist{i} j.*ones(1,round(TP_normalized(j,i).*10000))];
    end
    HLSDist{i}(end+1:10000)=2;
end

%generate LLS distributions
LLSDist1 = cell(10,1);LLSDist2 = cell(10,5);
for i = 1:10
    for j = 1:5
        for k = 1:5
            LLSDist2{i,j} = [LLSDist2{i,j} j.*ones(1,round(TP_LL_normalized{i}(k,j).*10000))];
        end
        LLSDist2{i,j}(end+1:10000)=2;
        LLSDist1{i} = [LLSDist1{i} j.*ones(1,round(LLSProbDist(i,j).*10000))];
    end
    LLSDist1{i}(end+1:10000)=2;
end

%generate tracks of HLS and LLS
state = zeros(length(initState),1000);state(:,1) = initState;
stateLL = zeros(length(initState),1000);
for i = 1:1000
    for HLS = 1:10
        currHLS = find(state(:,i)==HLS);
        nextHLSNdx = ceil(10000.*rand(length(currHLS),1));
        state(currHLS,i+1) = HLSDist{HLS}(nextHLSNdx);
        
        if i == 1
            curLLSNdx = ceil(10000.*rand(length(currHLS),1));
            stateLL(currHLS,1) = LLSDist1{HLS}(curLLSNdx);
        end
        
        sameState = currHLS(state(currHLS,i+1) == state(currHLS,i));
        diffState = setdiff(currHLS,sameState);
        if ~isempty(diffState)
            curLLSNdx = ceil(10000.*rand(length(diffState),1));
            for j = 1:length(diffState)
                stateLL(diffState(j),i+1) = LLSDist1{state(diffState(j),i+1)}(curLLSNdx(j));
            end
        end
        
        if ~isempty(sameState)
            for LLS = 1:5
                currLLS = find(stateLL(sameState,i) == LLS);
                nextLLSNdx = ceil(10000.*rand(length(currLLS),1));
                stateLL(sameState(currLLS),i+1) = LLSDist2{HLS,LLS}(nextLLSNdx);
            end
        end
    end
end
end

function [spdMat,curvMat] = synSpdCurvGen(state,stateLL,empLLSC)

%sample vPar and vPerp from LLS distributions
spdMat = zeros(size(state));curvMat = zeros(size(state));
for HLS = 1:10
    for LLS = 1:5
        tempHLSNdx = find(state == HLS);
        tempLLSNdx = tempHLSNdx(stateLL(tempHLSNdx) == LLS);
        currSCDist = empLLSC{HLS,LLS};
        SCNdx = ceil(rand(1,length(tempLLSNdx)).*size(currSCDist,2));
        curvMat(tempLLSNdx) = currSCDist(2,SCNdx);
        spdMat(tempLLSNdx) = currSCDist(1,SCNdx);
    end
end

end

function [spdMat,curvMat] = synSCGenRTStyle(state,stateLL,empLLSC,cond)

%sample vPar and vPerp from LLS distributions
spdMat = zeros(size(state));curvMat = zeros(size(state));
for HLS = 1:10
    for LLS = 1:5
        tempHLSNdx = find(state == HLS);
        tempLLSNdx = tempHLSNdx(stateLL(tempHLSNdx) == LLS);
        currSCDist = empLLSC{HLS,LLS};
        currAvgSC = mean(currSCDist,2);
        SCNdx = ceil(rand(1,length(tempLLSNdx)).*size(currSCDist,2));
        if strcmpi(cond,'constSpeed')
            spdMat(tempLLSNdx) = currAvgSC(1);
            curvMat(tempLLSNdx) = currSCDist(2,SCNdx);
        elseif strcmpi(cond,'constCurv')
            spdMat(tempLLSNdx) = currSCDist(1,SCNdx);
            curvMat(tempLLSNdx) = currAvgSC(2);
        end
    end
end

if strcmpi(cond,'constCurv')
    for i = 1:size(state,1)
        track = curvMat(i,:);
        transitionNdx = find(diff(track)~=0);
        startNdx = [2 transitionNdx+2];
        endNdx = [transitionNdx length(track)];
        for j = 1:length(startNdx)
            track(startNdx(j):endNdx(j))=0;
        end
        curvMat(i,:) = track;
    end
end

end

function [spdMat,curvMat] = vParPerp2SC(vParMat,vPerpMat)

%calculate speed and curvature from vPar and vPerp
spdMat = sqrt(vParMat.^2+vPerpMat.^2);
curvMat = atan(vPerpMat./vParMat);
curvature = -curvMat;

%the following portion mimicks parts of the preprocessing code
idx3=find(spdMat<0.03);
curvature(idx3)=0;
% removing curvature > 6. These are most likely to be some form of 0
% being 2pi
idx4=find(abs(curvature)>4);
curvature(idx4)=0;

idx5=find(abs(curvature)>2);
for pp=1:length(idx5)
    curvature(idx5(pp))=curvature(idx5(pp)-1);
end
for i = 1:size(curvature,1)
    curvature(i,:)=smooth(curvature(i,:),10);
end
curvMat = curvature;

end

function [xPosAll,yPosAll,distAll] = XYPosCalc(xTemp,yTemp,spdMat,curvMat)

xPosAll = zeros(size(spdMat));
yPosAll = zeros(size(spdMat));
xPosAll(:,1) = xTemp;
yPosAll(:,1) = yTemp;
currAngle = zeros(size(xPosAll,1),1);

for i = 1:1000
    currAngle = currAngle+curvMat(:,i);
    xPosAll(:,i+1) = xPosAll(:,i)+(spdMat(:,i).*cos(currAngle)./30)./3.2;
    yPosAll(:,i+1) = yPosAll(:,i)+(spdMat(:,i).*sin(currAngle)./30)./3.2;
    rPos = sqrt(xPosAll(:,i+1).^2+yPosAll(:,i+1).^2);
    
    ndx=find(rPos>1);
    xPosAll(ndx,i+1) = xPosAll(ndx,i);
    yPosAll(ndx,i+1) = yPosAll(ndx,i);
end
distAll = sqrt((xPosAll-xPosAll(:,1)).^2+(yPosAll-yPosAll(:,1)).^2);

end

% Plotting Functions
function HHMMPlottingFunction(spdMat,curvMat,allState,state,printFig)

%plot autocorrelation
figure(100)
acf = autocorr(spdMat(1,:),500);
subplot(2,1,1);hold on;plot(acf,'b')
acf = autocorr(curvMat(1,:),500);
subplot(2,1,2);hold on;plot(acf,'b')
legend({'Emp','HHMM'})
set(gcf,'Position',[849 49 824 918])
if printFig == 1
    print('-dpsc2',[['PDF_Files/ModelComparison'] '.ps'],'-loose','-append');
end

%plot HLS probabilities
HHMMHLSDist = zeros(1,10);empHLSDist = zeros(1,10);
HHMMHLSDistInd = zeros(1,10);empHLSDistInd = zeros(1,10);
for HLS = 1:10
    HHMMHLSDist(HLS) = sum(sum(allState==HLS));
    empHLSDist(HLS) = sum(sum(state==HLS));
    HHMMHLSDistInd(HLS) = sum(sum(allState(1,:)==HLS));
    empHLSDistInd(HLS) = sum(sum(state(1,:)==HLS));
end
HHMMHLSDist = HHMMHLSDist./sum(HHMMHLSDist);
empHLSDist = empHLSDist./sum(empHLSDist);
HHMMHLSDistInd = HHMMHLSDistInd./sum(HHMMHLSDistInd);
empHLSDistInd = empHLSDistInd./sum(empHLSDistInd);
figure(101)
subplot(2,1,1)
plot(empHLSDist,'r');hold on;plot(HHMMHLSDist,'b');legend({'Emp','HHMM'})
xlabel('HLS');ylabel('Probability');title('HLS Probability');ylim([0 0.3])
subplot(2,1,2)
plot(empHLSDistInd,'r');hold on;plot(HHMMHLSDistInd,'b');legend({'Emp','HHMM'})
xlabel('HLS');ylabel('Probability');title('HLS Probability');ylim([0 0.6])
set(gcf,'Position',[849 49 824 918])
if printFig == 1
    print('-dpsc2',[['PDF_Files/ModelComparison'] '.ps'],'-loose','-append');
end

end

function [] = plotSampTracks(T1,T2,T3,T4,c,printFig)

CF = 6;
% plot empirical and synthetic tracks
for fig = 1:2
    figure(fig)
    for i = 1:6
        subplot(3,2,i)
        hold on
        viscircles([0,0],1,'Color',[0.5 0.5 0.5])  
        %viscircles([0,0],1.5/3.2,'Color',[0.5 0.5 0.5],'LineStyle','--')

        plot(T1.x((fig-1).*CF+i,:),T1.y((fig-1).*CF+i,:),'Color',c(1,:));
        plot(T1.x((fig-1).*CF+i,[3,30,90,300]),T1.y((fig-1).*CF+i,[3,30,90,300]),'k.','MarkerSize',12);
        plot(T2.x((fig-1).*CF+i,:),T2.y((fig-1).*CF+i,:),'Color',c(2,:));
        plot(T3.x((fig-1).*CF+i,:),T3.y((fig-1).*CF+i,:),'Color',c(3,:));
%         plot(T4.x((fig-1).*CF+i,:),T4.y((fig-1).*CF+i,:),'Color',c(4,:));
        
        plot(T1.x((fig-1).*CF+i,1),T1.y((fig-1).*CF+i,1),'k*')
        
        axis([-1 1 -1 1])
        title(num2str((fig-1).*CF+i))
    end
    set(gcf,'Position',[849 49 824 918])
    if printFig == 1
        print('-dpsc2',[['PDF_Files/ModelComparison'] '.ps'],'-loose','-append');
    end
end

end

function [x,y,d,d2] = plotModelComparisons(distAll1,distAll2,distAll3,distAll4,color)
dist2Consider = [1,2,3,15,30,90,100,300,900];
%dist2Consider = [1,2,3,15:15:900];
x = cell(1,length(dist2Consider));y = cell(3,length(dist2Consider));
d = zeros(3,length(dist2Consider));d2 = zeros(3,length(dist2Consider));
for i = 1:length(dist2Consider)
    if i<10
        subplot(3,3,i)
    end
    hold on
    curr1 = distAll1(:,dist2Consider(i)+1);
    curr2 = distAll2(:,dist2Consider(i)+1);
    curr3 = distAll3(:,dist2Consider(i)+1);
    curr4 = distAll4(:,dist2Consider(i)+1);
    
    x_values = linspace(0,max([curr1;curr2;curr3;curr4]), 50);
    pd = fitdist(curr1,'Kernel','Kernel','epanechnikov');
    y_1 = pdf(pd,x_values);
    pd = fitdist(curr2,'Kernel','Kernel','epanechnikov');
    y_2 = pdf(pd,x_values);
    pd = fitdist(curr3,'Kernel','Kernel','epanechnikov');
    y_3 = pdf(pd,x_values);
    pd = fitdist(curr4,'Kernel','Kernel','epanechnikov');
    y_4 = pdf(pd,x_values);
    
    plot(x_values,y_1,'Color',color(1,:))
    plot(x_values,y_2,'Color',color(2,:))
    plot(x_values,y_3,'Color',color(3,:))
    plot(x_values,y_4,'Color',color(4,:))
    title(['t = ' num2str(dist2Consider(i))])
    y_1Norm = y_1'./sum(y_1);
    temp = (y_1Norm == 0);
    y_1Norm(temp) = y_1Norm(temp)+eps;
    y_1Norm(~temp) = y_1Norm(~temp)-eps.*sum(temp)./(sum(~temp));
    
    y_2Norm = y_2'./sum(y_2);
    temp = (y_2Norm == 0);
    y_2Norm(temp) = y_2Norm(temp)+eps;
    y_2Norm(~temp) = y_2Norm(~temp)-eps.*sum(temp)./(sum(~temp));
    
    y_3Norm = y_3'./sum(y_3);
    temp = (y_3Norm == 0);
    y_3Norm(temp) = y_3Norm(temp)+eps;
    y_3Norm(~temp) = y_3Norm(~temp)-eps.*sum(temp)./(sum(~temp));
    
    y_4Norm = y_4'./sum(y_4);
    temp = (y_4Norm == 0);
    y_4Norm(temp) = y_4Norm(temp)+eps;
    y_4Norm(~temp) = y_4Norm(~temp)-eps.*sum(temp)./(sum(~temp));
    
    d(1,i)=jensen_shannon_divergence(y_1Norm,y_2Norm);
    d(2,i)=jensen_shannon_divergence(y_1Norm,y_3Norm);
    d(3,i)=jensen_shannon_divergence(y_1Norm,y_4Norm);

    d2(1,i)=sum(sqrt(y_1Norm.^2+y_2Norm.^2));
    d2(2,i)=sum(sqrt(y_1Norm.^2+y_3Norm.^2));
    d2(3,i)=sum(sqrt(y_1Norm.^2+y_4Norm.^2));
    
    x{i} = x_values;
    y{1,i} = y_1;
    y{2,i} = y_2;
    y{4,i} = y_3;
    y{3,i} = y_4;
end

end

