function Leader=Leader_Selection(Archive,beta)

    GridIndices=[Archive.GridIndex];
    
    occ_cell_index=unique(GridIndices);
    
    occ_cell_member_count=zeros(size(occ_cell_index));

    m=numel(occ_cell_index);
    for k=1:m
        occ_cell_member_count(k)=sum(GridIndices==occ_cell_index(k));
    end
    
    p=occ_cell_member_count.^(-beta);
    p=p/sum(p);
    
    selected_cell_index=occ_cell_index(find((rand)<=cumsum(p),1,'first'));
    
    selected_cell_members=find(GridIndices==selected_cell_index);
    
    selected_member_index=randi([1 numel(selected_cell_members)]);
    
    Leader=Archive(selected_cell_members(selected_member_index));
end
