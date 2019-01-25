function [data_time,model,S,likely_high_state_sorted,HHMM_by_State_sorted,...
    HHMM_lowlevel_sorted,likely_high_state_by_fly,clusters_to_consider,ndx,ndxLL,dur_sorted] = Processing(options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% The Following Section is for loading and preprocessing %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code loads the data and model file(s) and apply corrections
% to when the fly is in the odor ring data(12,:) and the stim on data(11,:)
%
% Files needed to run code:
%       1.) .mat file for data
%       2.) post_hoc_preprocess_LT.m (located under '\Support Files\HHMM_Support')
%       3.) HHMM and HMM model files (located under '\Support Files\HHMM_Support')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(options.datFile,'model','data')
fig_count=0;
% Update locations where fly is in odor ring
data_old = data;
data = post_hoc_preprocess_LT(data_old);
% add in time step to row 17 of data
data_time = data;
for i = 1:size(data_time,2)
    data_time{i} = [data{i}; 1:1:size(data{i},2)];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% The Following Sections are for processing %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code has 3 sections.
%       Section 1: Parsing data into HLS and store them as HLS and LLS tracks
%           Main Outputs: Cell arrays HHMM_by_State, HHMM_lowlevel
%       Section 2: Sorting and saving speed curvature for all HLS
%           Main Outputs: Structure S.HL
%       Section 3: Sorting HL states by speed and curvature
%           Main Outputs: Array likely_high_state_sorted
%       Section 4: Sorting and saving speed curvature for LLS
%           Main Outputs: Structure S.LL, S.data
%
% Outputs:
%       HHMM_by_State: Kinematics for each HLS track
%       HHMM_lowlevel: LLS compositions for each HLS track
%       likely_high_state_sorted: HLS for every time instance
%       S: Structure detailing all characteristics of each track
%           S.data: kinematics
%           S.HL: speed and curvature for all HLS (z, no z, obs)
%           S.LL: .avgspd, .avgcurv, .allTrack
%                   .allTrack: Highest likely LLS trajectory for HLS track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% parsing data into high level states. Each state is stored as a cell
number_of_flies_in_a_cluster=sum(model.p,2);                                % essentially model.p tells which model for which fly
clusters_to_consider=find(number_of_flies_in_a_cluster>=10);                  % find models that fit at least  flies
HHMM_by_State = cell(length(clusters_to_consider),sum(model.NA),model.Qdim,1);
HHMM_lowlevel = cell(length(clusters_to_consider),sum(model.NA),model.Qdim,1);
for i=1:length(clusters_to_consider)
    cluster_number=clusters_to_consider(i);
    HHMMcluster=model.HHMMs{1,cluster_number}.Qp;
    fly_no=find(model.p(cluster_number,:)>0.5);                             % this finds all the flies which are fit by a given model
    durations = cell(model.Qdim,1);
    for ii=1:length(fly_no)
        HHMMforthefly=HHMMcluster{1,fly_no(ii)};
        [prob,likely_high_state]=max(HHMMforthefly);
        conditional=find(prob<0.85);                                        % includes only points with >85% probability
        likely_high_state(conditional)=0;                                   % changes the states with p<0.85 to 0
        tsig = likely_high_state > 0;
        
        %This looks at places where the state changes to another state for
        %less than 5 steps before changing back to initial state
        dsig = diff([1 tsig 1]);
        startIndex = find(dsig < 0);
        endIndex = find(dsig > 0)-1;
        dur = endIndex-startIndex+1;
        stringIndex = (dur <= 5);
        startIndex = startIndex(stringIndex);
        endIndex = endIndex(stringIndex);
        
        indices = zeros(1,max(endIndex)+1);
        indices(startIndex) = 1;
        indices(endIndex+1) = indices(endIndex+1)-1;
        indices = find(cumsum(indices));
        
        for ij=1:length(indices)
            if indices(ij)>1
                likely_high_state(indices(ij))=likely_high_state(indices(ij)-1);
            end
        end
        
        %This removes regions where a state only occurs for one instance
        %and changes it to the state at the next time point
        find(diff(likely_high_state)==1);
        a = find(~(diff([0, likely_high_state])==0));
        aa = find(diff(a)==1);
        for jj=length(aa):-1:1
            likely_high_state(a(aa(jj)))=likely_high_state(a(aa(jj))+1);
        end
        
        likely_high_state_by_fly(fly_no(ii),1:size(data_time{1,1},2))=likely_high_state;
        statetransitions=find(diff(likely_high_state)~=0);                  % calculates the time at which states transitions occur
        state_begin=unique([1 statetransitions(1:end-1)+1]);
        state_end=unique(statetransitions);
        counter=zeros(model.Qdim,1);
        for jj=1:length(statetransitions)                                   % iterating through each high-level state
            index=fly_no(ii);                                               % fly number
            state_index=likely_high_state(state_begin(jj));
            if state_index==0
            else
                durations{state_index} = [durations{state_index} state_end(jj)-state_begin(jj)+1];
                
                counter(state_index)=counter(state_index)+1;
                if state_begin(jj)>10
                    kinematics=data_time{1,index}(:,(state_begin(jj)-10):state_end(jj));
                    llbegin=(state_index-1)*model.dim+1;
                    llend=state_index*model.dim;
                    lowlevelstate=model.HHMMs{1,1}.p{1,index}(llbegin:llend,(state_begin(jj)-10):state_end(jj));
                    kinematics=num2cell(kinematics,[1,2]);
                    lowlevelstate=num2cell(lowlevelstate,[1,2]);
                    HHMM_by_State(i,ii,state_index,counter(state_index))=kinematics;
                    HHMM_lowlevel(i,ii,state_index,counter(state_index))=lowlevelstate;
                end
            end
        end
    end
end

%% Sorting and saving speed curvature for all high level states
if exist('fig_count','var')==0
    fig_count=1;
else
    fig_count=fig_count+1;
end

%Initialize matrices
S.HL = [];
ratio = zeros(size(HHMM_by_State,1),model.Qdim);

avgspeed = cell(size(HHMM_by_State,1),model.Qdim);
avgcurvature = cell(size(HHMM_by_State,1),model.Qdim);
speed_zT = cell(size(HHMM_by_State,1),model.Qdim);                          %z-scored
curvature_zT = cell(size(HHMM_by_State,1),model.Qdim);                      %z-scored
v_par_T = cell(size(HHMM_by_State,1),model.Qdim);
v_perp_T = cell(size(HHMM_by_State,1),model.Qdim);

S.HL.avgspd = cell(size(HHMM_by_State,1),model.Qdim);
S.HL.avgcurv = cell(size(HHMM_by_State,1),model.Qdim);
S.HL.spd_z = cell(size(HHMM_by_State,1),model.Qdim);                        %z-scored
S.HL.curv_z = cell(size(HHMM_by_State,1),model.Qdim);                       %z-scored
S.HL.v_par = cell(size(HHMM_by_State,1),model.Qdim);
S.HL.v_perp = cell(size(HHMM_by_State,1),model.Qdim);

ndx = cell(1,size(HHMM_by_State,1));

%Loop through each cluster and state
for k=1:size(HHMM_by_State,1)
    currentcluster=squeeze(HHMM_by_State(k,:,:,:));
    for kk=1:size(currentcluster,2)
        currentstate=squeeze(currentcluster(:,kk,:));                       % currentstate = all instances of state kk in cluster k
        speed=[];
        meanspeed=[];
        curvature=[];
        meancurvature=[];
        speed_z = [];
        meanspeed_z = [];
        curvature_z = [];
        meancurvature_z = [];
        v_perp = [];
        mean_v_perp = [];
        v_par = [];
        mean_v_par = [];
        
        currentstatenew = currentstate(~cellfun(@isempty, currentstate));
        if isempty(currentstatenew)==1
        else
            [a,duration]=cellfun(@size, currentstatenew);
            % Only include tracks of HLS longer than 25 instances
            longduration=find(duration>25);                                 
            number=length(longduration);
            for kkk=1:number
                index=longduration(kkk);
                currenttrack=currentstatenew{index};
                
                speed=[speed currenttrack(7,:)];
                meanspeed=[meanspeed nanmean(currenttrack(7,:))];
                
                curvature=[curvature currenttrack(8,:)];
                meancurvature=[meancurvature nanmean(currenttrack(8,:))];
                
                speed_z = [speed_z currenttrack(13,:)];
                meanspeed_z = [meanspeed_z nanmean(currenttrack(13,:))];
                
                curvature_z = [curvature_z currenttrack(14,:)];
                meancurvature_z=[meancurvature_z nanmean(currenttrack(14,:))];
                
                v_par = [v_par currenttrack(15,11:end)];
                mean_v_par=[mean_v_par nanmean(currenttrack(15,11:end))];
                
                v_perp = [v_perp currenttrack(16,11:end)];
                mean_v_perp = [mean_v_perp nanmean(currenttrack(16,11:end))];
                
            end
            avgspeed{k,kk}=meanspeed;
            avgcurvature{k,kk}=meancurvature;
            speed_zT{k,kk}=meanspeed_z;
            curvature_zT{k,kk}=meancurvature_z;
            v_par_T{k,kk} = mean_v_par;
            v_perp_T{k,kk} = mean_v_perp;
            
            % Define the ratio as speed over curvature
            ratio(k,kk) = mean(meanspeed)./std(meancurvature);              
        end
    end
    
    %Sort the ratio from low to high and store to another matrix called ndx
    ratio_temp = ratio(k,1:size(HHMM_by_State,3));
    [~,ndxTemp] = sort(ratio_temp);
    ndx{k}=ndxTemp;
    
    %Sort the HL and LL states based on the ratio ranking in the current
    %cluster
    HHMM_by_State_sorted(k,:,:,:) = HHMM_by_State(k,:,ndx{k}(1:size(HHMM_by_State,3)),:);
    HHMM_lowlevel_sorted(k,:,:,:) = HHMM_lowlevel(k,:,ndx{k}(1:size(HHMM_by_State,3)),:);
    
    %Sort the avgspeed and avgcurvature based on the ndx
    for kk = 1:size(currentcluster,2)
        S.HL.avgspd{k,kk}=avgspeed{k,ndxTemp(kk)};
        S.HL.avgcurv{k,kk}=avgcurvature{k,ndxTemp(kk)};
        
        S.HL.spd_z{k,kk}=speed_zT{k,ndxTemp(kk)};
        S.HL.curv_z{k,kk}=curvature_zT{k,ndxTemp(kk)};
        
        S.HL.v_par{k,kk}=v_par_T{k,ndxTemp(kk)};
        S.HL.v_perp{k,kk}=v_perp_T{k,ndxTemp(kk)};
    end
end

%% Sorting HL states by speed and curvature
likely_high_state_sorted = likely_high_state_by_fly;
for i=1:length(clusters_to_consider)                                        % looping through the clusters
    cluster=clusters_to_consider(i);
    %fly_no=find(model.p(cluster,:)>0.5);
    fig_count=fig_count+1;                                                       % first entry of the fly after odor on
    fly_no=find(model.p(cluster,:)>0.5);
    dur_sorted = cell(model.Qdim,1);
    for ii=1:length(fly_no)                                                 % within that cluster finding the high-level states
        high_probs=likely_high_state_by_fly(fly_no(ii),:);
        for iii = 1:length(ndx{i})
            likely_high_state_sorted(fly_no(ii),high_probs == ndx{i}(iii)) =iii;% sort the HL states based on ndx
        end
    end
    for iii = 1:length(ndx{i})
        dur_sorted(iii) = durations(ndx{i}(iii));
    end
end

%% Sorting and saving speed curvature for low level states
S.LL = [];
S.LL.avgspd.allTrack = cell(size(HHMM_by_State,1),model.Qdim,model.dim);
S.LL.avgcurv.allTrack = cell(size(HHMM_by_State,1),model.Qdim,model.dim);
S.LL.avgspd.noNaN = cell(size(HHMM_by_State,1),model.Qdim,model.dim);
S.LL.avgcurv.noNaN = cell(size(HHMM_by_State,1),model.Qdim,model.dim);
S.LL.allTrack = cell(size(HHMM_by_State,1));
ndxLL = cell(size(HHMM_by_State_sorted,1), size(HHMM_by_State_sorted,3));

for k=1:size(HHMM_by_State_sorted,1)
    currentcluster=squeeze(HHMM_by_State_sorted(k,:,:,:));
    currentlowlevelcluster=squeeze(HHMM_lowlevel_sorted(k,:,:,:));          % access sorted LL states
    for kk=1:size(HHMM_by_State_sorted,3)
        speed={};
        meanspeed={};
        curvature={};
        meancurvature={};
        currentstate=squeeze(currentcluster(:,kk,:));
        currentstatenew = currentstate(~cellfun(@isempty, currentstate));
        
        if isempty(currentstatenew)==1
        else
            lowlevelstate=squeeze(currentlowlevelcluster(:,kk,:));
            [a,duration]=cellfun(@size, currentstatenew);
            longduration=find(duration>25);                                 % Only include tracks longer than 20 instances
            lowlevelstatenew = lowlevelstate(~cellfun(@isempty,lowlevelstate));
            [a,duration]=cellfun(@size, currentstatenew);
            longduration=find(duration>25);
            number=length(longduration);
            for kkk=1:number                                                % loop through all tracks in a HL state
                index=longduration(kkk);
                currenttrack=currentstatenew{index};
                currentlowlevelstate=lowlevelstatenew{index};
                [~,currentlowlevel]=max(currentlowlevelstate);
                S.LL.allTrack{k}{kk}{kkk} = currentlowlevel;
                S.data{k}{kk}{kkk} = currenttrack;
                for kkkk=1:model.dim                                        % loop through all LL states in the HL state
                    idx=find(currentlowlevel==kkkk);                        % idx all tracts of LL state kkkk in HL state kk
                    speed{kkk,kkkk}={currenttrack(13,idx)};
                    meanspeed{kkk,kkkk}= {nanmean(currenttrack(13,idx))};
                    curvature{kkk,kkkk}={ currenttrack(14,idx)};
                    meancurvature{kkk,kkkk}={nanmean(currenttrack(14,idx))};
                end
            end
            counts_tot = zeros(20);
            for kkkk=1:model.dim
                speed_all=[];
                meanspeed_all=[];
                curvature_all=[];
                meancurvature_all=[];
                for kkk=1:number
                    a=cell2mat(speed{kkk,kkkk});
                    b=cell2mat(meanspeed{kkk,kkkk});
                    c=cell2mat(curvature{kkk,kkkk});
                    d=cell2mat(meancurvature{kkk,kkkk});
                    speed_all=[speed_all a];
                    meanspeed_all=[meanspeed_all b];                        % putting mean speed into a array rather than a cell array
                    curvature_all=[curvature_all c];
                    meancurvature_all=[meancurvature_all d];                % putting mean curvature into array rather than cell array
                end
                avgspeedLL{k,kk,kkkk}=meanspeed_all;
                avgcurvatureLL{k,kk,kkkk}=meancurvature_all;
                meanspeed_all = meanspeed_all(~isnan(meanspeed_all));
                meancurvature_all = meancurvature_all(~isnan(meancurvature_all));
                ratioLL(k,kk,kkkk) = mean(meanspeed_all)./std(meancurvature_all); % Define the ratio as speed over curvature
            end
            
            %Sort the LL ratio from low to high and store to another matrix
            %called ndxLL
            ratio_temp = ratioLL(k,kk,:);
            [~,ndxTemp] = sort(ratio_temp);
            ndxLL{k,kk} = ndxTemp;
            
            %Sort the LL avgspeed and avgcurvature based on the ndxLL
            for kkkk = 1:model.dim
                tempspdSorted = avgspeedLL{k,kk,ndxLL{k,kk}(kkkk)};
                tempcurvSorted = avgcurvatureLL{k,kk,ndxLL{k,kk}(kkkk)};
                
                S.LL.avgspd.allTrack{k,kk,kkkk}=tempspdSorted;
                S.LL.avgcurv.allTrack{k,kk,kkkk}=tempcurvSorted;
                S.LL.avgspd.noNaN{k,kk,kkkk}=tempspdSorted...
                    (~isnan(tempspdSorted));
                S.LL.avgcurv.noNaN{k,kk,kkkk}=tempspdSorted...
                    (~isnan(tempspdSorted));
            end
        end
    end
end

end