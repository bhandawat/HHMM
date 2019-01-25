function [] = plotXYPositionTracks(data,model,likely_high_state_sorted,clusters_to_consider,fig_title)
if isempty(clusters_to_consider)
    number_of_flies_in_a_cluster=sum(model.p,2);                            % essentially model.p tells which model for which fly
    clusters_to_consider=find(number_of_flies_in_a_cluster>=10);            % find models that fit at least  flies
end

colorMap=varycolor(model.Qdim+1);
cc=cell(1,length(clusters_to_consider));
for i = 1:length(clusters_to_consider)
    cc{i}=colorMap;
    cc{i}(end,:) = [1 1 1];
end
nFlyPerPage = 6;

for i=1:length(clusters_to_consider)                                        % looping through the clusters
    cluster=clusters_to_consider(i);
    first_entry1=[];                                                        % first entry of the fly after odor on
    fly_no = find(model.p(cluster,:)>0.75);
    for ii=1:length(fly_no)                                                 % within that cluster finding the high-level states
        rowNdx = mod(ii,nFlyPerPage);
        rowNdx(rowNdx==0) = nFlyPerPage;
        figNdx = ceil(ii./nFlyPerPage);
        
        likely_high_state = likely_high_state_sorted(fly_no(ii),:);
        likely_high_state(likely_high_state == 0) = 11;
        
        d1 = data{fly_no(ii)};
        
        idx{1}=find(d1(11,1:end-1)<1 & d1(12,1:end-1)<1);        % odor off and fly outside
        idx{2}=find(d1(11,1:end-1)<1 & d1(12,1:end-1)==1);       % odor off and fly inside
        idx{3}=find(d1(11,1:end-1)==1 & d1(12,1:end-1)<1);        % odor on and fly outside
        idx{4}=find(d1(11,1:end-1)==1 & d1(12,1:end-1)==1);       % odor on and fly inside
        titles = {'b_o','b_i','d_o','d_i'};
        
        figure(figNdx)
        for iii = 1:4
            subplot(nFlyPerPage,4,(rowNdx-1)*4+iii)
            scatter(d1(1,idx{iii}),d1(2,idx{iii}),1,[cc{i}(likely_high_state(idx{iii}),:)]);hold on
            viscircles([0,0],1.5/3.2,'LineStyle','--','EdgeColor',[0.7 0.7 0.7]);
            viscircles([0,0],1,'EdgeColor',[0.7 0.7 0.7]);
            xlim([-1 1]);ylim([-1 1]);
            if rowNdx == 1
                title(titles{iii})
            end
            if iii ==1
                ylabel(['Fly ' num2str(fly_no(ii))])
            end
        end
    end
    for ii = 1:figNdx
        figure(ii);set(gcf,'position',[849 49 824 918])
        print('-dpsc2',[fig_title '.ps'],'-loose','-append');
    end
    close all
end

end