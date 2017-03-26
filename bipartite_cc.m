function [cell_sub_nets]=bipartite_cc(map)

%********************
% map is a binary bipartite graph
[m,n]=size(map);

%********************
% making a network with m+n nodes
% first group 1:m, and second group m+1:m+n

map2=sparse(m+n,m+n,0,m+n,m+n);
[x,y]=find(map~=0);
for i=1:length(x)
    map2(x(i),y(i)+m)=1;
    map2(y(i)+m,x(i))=1;
end

%********************
% computing connected components
[labels rts] = graph_connected_components(map2);
len=length(unique(labels));

%********************
% computing node labels for each connected component
cell_sub_nets={};
flag=1;
for i=1:len
    ind=find(labels==i);
    if length(ind)>=2
        ind_ind=find(ind<=m);
        row_nodes=ind(ind_ind);
        col_nodes=setdiff(ind,row_nodes)-m;
        cell_sub_nets{flag,1}=[row_nodes];
        cell_sub_nets{flag,2}=[col_nodes];
        flag=flag+1;
    end
end


