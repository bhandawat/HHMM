function [model,data,HMMsortNdx,likely_state_by_fly_sorted,dur] = ProcessingHMM(options)
%clear
%load('Data/Sept9_WTACV0_HMM_50_var2.mat')
%load('Data/Sept10_WTACV0_HMM_50_var2.mat')
load(options.datFileHMM)
nFlys = 34;
likely_state_by_fly = zeros(nFlys,10799);
dur = cell(model.dim,1);
Track = cell(model.dim,nFlys);

for fly = 1:nFlys
    d1 = data{fly};
    [p,ndx] = max(model.HMMs{1}.p{fly});  
    
%-------------------------------------------------------------------------
    conditional=find(p<0.85);                             % includes only points with >85% probability
    ndx(conditional)=0;                                   % changes the states with p<0.85 to 0
    tsig = ndx > 0;                                       %# Using eps as the threshold
    
    %This looks at places where the state changes to another state for
    %less than 5 steps before changing back to initial state
    dsig = diff([1 tsig 1]);
    startIndex = find(dsig < 0);
    endIndex = find(dsig > 0)-1;
    duration = endIndex-startIndex+1;
    stringIndex = (duration <= 5);
    startIndex = startIndex(stringIndex);
    endIndex = endIndex(stringIndex);
    
    indices = zeros(1,max(endIndex)+1);
    indices(startIndex) = 1;
    indices(endIndex+1) = indices(endIndex+1)-1;
    indices = find(cumsum(indices));
    
    for ij=1:length(indices)
        if indices(ij)>1
            ndx(indices(ij))=ndx(indices(ij)-1);
        end
    end
    
    %This removes regions where a state only occurs for one instance
    %and changes it to the state at the next time point
    find(diff(ndx)==1);
    a = find(~(diff([0, ndx])==0));
    aa = find(diff(a)==1);
    for jj=length(aa):-1:1
        ndx(a(aa(jj)))=ndx(a(aa(jj))+1);
    end
%-------------------------------------------------------------------------
    likely_state_by_fly(fly,:) = ndx;

    
    for state = 1:model.dim
        temp = likely_state_by_fly(fly,:)==state;
        
        startNdx = find(diff([false temp])==1);
        endNdx = find(diff([temp false])==-1);
        
        if length(startNdx)>length(endNdx)
            endNdx = [endNdx 10799];
        end
        
        if length(endNdx)>length(startNdx)
            startNdx = [1 startNdx];
        end
            
        dur{state} = [dur{state} endNdx-startNdx+1];
        for t = 1:length(startNdx)
            Track{state,fly}{t} = d1(:,startNdx(t):endNdx(t));
            if endNdx(t)-startNdx(t)>-1
                longTrack{state,fly}{t} = d1(:,startNdx(t):endNdx(t)); 
            end
            
        end
    end
end

avgvPar = ones(1,model.dim);stdvPerp = zeros(1,model.dim)+eps;
for state=1:model.dim
    allTrack = longTrack(state,:);
    allTrack = cat(2, allTrack{:});
    vPar = [];vPerp = [];
    if ~isempty(allTrack)
        allTrack = allTrack(~cellfun('isempty',allTrack));
        
        for t = 1:length(allTrack)
            vPar = [vPar nanmean(allTrack{t}(7,:))];
            vPerp = [vPerp nanmean(allTrack{t}(8,:))];
        end
        
        avgvPar(state)=mean(vPar);
        stdvPerp(state)=std(vPerp);
    end
end

ratio = avgvPar./stdvPerp;

[HMMsortRatio,HMMsortNdx] = sort(ratio);

likely_state_by_fly_sorted = likely_state_by_fly;                                                     % first entry of the fly after odor on
fly_no = 1:nFlys;
for ii=1:length(fly_no)                                                 % within that cluster finding the high-level states
    high_probs=likely_state_by_fly(fly_no(ii),:);
    for iii = 1:length(HMMsortNdx)
        likely_state_by_fly_sorted(fly_no(ii),high_probs == HMMsortNdx(iii)) =iii;% sort the HL states based on ndx
    end
end

% calculate duration
durHist = zeros(model.dim,2000);
maxTrack = zeros(model.dim,1);
for state = 1:model.dim
    temp = dur{state};
    if ~isempty(temp)
        n = histcounts(temp,max(temp));
        durHist(state,1:max(temp)) = n;
        maxTrack(state) = max(temp);
    end
end

% load('tempDurHist.mat');
% 
% figure;
% subplot(2,1,1)
% plot([1:1:length(durHist)]/30,sum(durHist)./sum(durHist(:)),'g')
% hold on;
% plot([1:1:length(HHMMDurHist)]/30,sum(HHMMDurHist)./sum(HHMMDurHist(:)),'k')
% xlabel('Duration (s)')
% ylabel('% tracks')
% xlim([0 1.5]);ylim([0 0.25])
% legend({'HMM','HHMM'})
% 
% subplot(2,1,2)
% plot([0:1:length(durHist)]/30,[1 1-cumsum(sum(durHist))./sum(durHist(:))],'g')
% hold on;
% plot([0:1:length(HHMMDurHist)]/30,[1 1-cumsum(sum(HHMMDurHist))./sum(HHMMDurHist(:))],'k')
% xlabel('Duration (s)')
% ylabel('1-cumsum(% tracks)')
% xlim([0 1.5]);ylim([0 1])
% legend({'HMM','HHMM'})
end


