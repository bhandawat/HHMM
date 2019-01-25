function [data_new, first_entry] = post_hoc_preprocess_LT(data)
data_new = data;
first_entry = zeros(1,length(data));
for ii = 1:length(data)
    x = data{ii}(1,:);
    y = data{ii}(2,:);
    R = sqrt(x.^2+y.^2);
    
    % odor on or off
    idx=find(R<(1.5/3.2));
    inside_ring = zeros(1,length(data{ii}));
    inside_ring(idx) = 1;
    
    data_new{ii}(12,:) = inside_ring;                                       %inside the odor ring
    %data_new{ii}(11,:) = [zeros(1,ceil(length(R)/2)), ones(1,ceil(length(R)/2)-1)]; % stim on
    
    idx = idx(find(idx>=length(R)/2));
    first_entry(ii) = idx(1);
    data_new{ii}(11,:) = [zeros(1,idx(1)-1), ones(1,ceil(length(R))-idx(1)+1)]; % stim on
    
    
end

end