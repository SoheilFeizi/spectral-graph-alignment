function map11_cell=EigenAlign(N,M,map,gamma_vec)


%--------------------------------------------------------------------------
% Spectral alignment of graphs
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% USAGE:
%    map11=EigenAlign(N,M,map,gamma)
%
% INPUT ARGUMENTS:
% N and M       Input networks, networks are binary square matrices,
%               possibly asymmetric.
%
% map:          provides a (many-to-many) binary mapping edges across graphs.
%               If two graphs have n1 and n2 nodes, map is a n1 by n2
%               binary matrix. Use map=ones(n1,n2) as default

% gamma_vec: a vector of input gamma parameters that controls the
% match/mismatch weight in the objective function of the network alignment
% optimization. Each element of the vector should be in [0,1/2)

% OUTPUT ARGUMENTS:
% map11_cell:provides a one-to-one mapping across networks for each gamma parameter.


% Beta version: 26-3-2017
% You can use it for non-commercial purposes
% All rights are reserved

% Please cite:
%  Spectral Alignment of Graphs
%  By: Soheil Feizi, Gerald Quon, Mariana Mendoza, Muriel Medard, Manolis Kellis, Ali Jadbabaie
%  arXiv: 1602.04181


if length(find(gamma_vec>=1/2 | gamma_vec<0))~=0
    disp('Error: gamma should be in [0,1/2)')
    return;
end

alpha_vec=(1./gamma_vec)-1;
eps=0.001;

%********************
% constructing alignment network
G=prep_align_net(full(N),full(M),full(map));
G=full(G);
%********************
ind_nz=find(sum(abs(G))~=0); % non-isolated nodes
map11_cell=cell(1,length(gamma_vec));
for gamma_iter=1:length(gamma_vec)
    disp(['running eigenalign:',num2str(gamma_iter),'-th gamma value'])
    alpha=alpha_vec(gamma_iter);
    G2=G(ind_nz,ind_nz);
    
    if alpha~=Inf
        G2(G2==1)=alpha+eps;
        G2(G2==-1)=eps;
        G2(G2==0)=1+eps;
    else
        % gamma=0, only matches
        G2(G2==1)=1;
        G2(G2==-1)=0;
        G2(G2==0)=0;
    end
    
    %disp('eigen decomposition')
    
    [u2,d2]=eigs(G2,1,'LA');
    u2=double(u2);
    
    %
    u2=abs(u2);
    u=zeros(size(G,1),1);
    u(ind_nz)=u2;
    
    [x_ind,y_ind]=find(map~=0);
    
    map_r=map*0;
    for i=1:length(x_ind)
        map_r(x_ind(i),y_ind(i))=u(i);
    end
    
    
    %********************
    % compute 1-1 mapping using max weight bipartite matching
    % to solve LP effeciently, we split the problem to connected components
    % over the bipartite graph and solve LP over each subproblem.
    % disp('LP')
    
    map11=map*0;
    
    cell_sub_nets=bipartite_cc(map);
    for i=1:size(cell_sub_nets,1)
        
        ind1=cell_sub_nets{i,1};
        ind2=cell_sub_nets{i,2};
        
        if length(ind1)>0 && length(ind2)>0
            if length(ind1)*length(ind2)==1
                map11(ind1,ind2)=1;
            else
                map_sub=map_r(ind1,ind2);
                %****
                map_sub2=-map_sub;
                map_sub2(map_sub==0)=inf;
                [map_sub_11,~] = Hungarian(map_sub2);
                map11(ind1,ind2)=(map_sub_11)>0;
            end
        end
    end
    
    map11_cell{gamma_iter}=map11;
end