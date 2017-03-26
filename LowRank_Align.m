function [map11_LA]=LowRank_Align(N1,N2,map,k_rank,gamma)

%--------------------------------------------------------------------------
% Spectral alignment of graphs
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% USAGE:
%    map11=LowRank_Align(N1,N2,map,kk,gamma)
%
% INPUT ARGUMENTS:
% N1 and N2:     Input graphs, graphss are binary, square and symmetric
% matrices.
%
% map:          provides a (many-to-many) binary mapping edges across graphs.
%               If two graphs have n1 and n2 nodes, map is a n1 by n2
%               binary matrix. Use map=ones(n1,n2) as default

% k_rank is the number of top eigenvectors to consider

% gamma: input gamma parameter that controls the
% match/mismatch weight in the objective function of the network alignment
% optimization. It should be in [0,1/2)

% OUTPUT ARGUMENTS:
% map11_cell:provides a one-to-one mapping across graphs

% Beta version: 26-3-2017
% You can use it for non-commercial purposes
% All rights are reserved

% Please cite:
%  Spectral Alignment of Graphs
%  By: Soheil Feizi, Gerald Quon, Mariana Mendoza, Muriel Medard, Manolis Kellis, Ali Jadbabaie
%  arXiv: 1602.04181


N1=double(N1);
N2=double(N2);

map=(map~=0)*1.0;

n1=size(N1,1);
n2=size(N2,1);

%*******************************
[U1,D1]=eigs(N1-gamma,k_rank,'LA');
[U2,D2]=eigs(N2-gamma,k_rank,'LA');


[~,I1]=sort(diag(D1),'descend');
U11=U1(:,I1);
D11=D1(I1,I1);

[~,I2]=sort(diag(D2),'descend');
U22=U2(:,I2);
D22=D2(I2,I2);

%***********************************
% making majority of eigenvector elements positive
% bc if we multiply -1, it wont change it

U11_p=U11;
U22_p=U22;

for j=1:k_rank
    if sum(U11(:,j)>=0)<(n1/2)
        U11_p(:,j)=-U11(:,j);
    end
    
    if sum(U22(:,j)>=0)<(n2/2)
        U22_p(:,j)=-U22(:,j);
    end
end

U11_p=double(U11_p);
U22_p=double(U22_p);

%*******************************
flip_mat = dec2bin(0:2^(k_rank)-1)-'0';
map11_cell=cell(1,size(flip_mat,1));

for ii=1:size(flip_mat,1)
    flip_vec=2*flip_mat(ii,:)-1;
    
    X_OLPM=zeros(n1,n2);
    for j=1:k_rank
        X_OLPM=X_OLPM+abs(D11(j,j)*D22(j,j))*flip_vec(j)*U11_p(:,j)*U22_p(:,j)';
        
    end
    X_OLPM(map==0)=-Inf;
    X_OLPM=double(X_OLPM);
    %********************
    % compute 1-1 mapping using max weight bipartite matching
    XX_OLPM=max_matching_Hungarian(X_OLPM);
    map11_cell{ii}=XX_OLPM;
end

%*******************************
% pick the one that leads to the highest objective
k_temp=size(map11_cell,2);
match_mismatch_temp=zeros(2,k_temp);
for i=1:k_temp
    [num_matches,num_mismatches]=eval_rate_match_mismatch2(N1,N2,map11_cell{i});
    match_mismatch_temp(:,i)=[num_matches;num_mismatches];
end

[~,ind_t]=max(match_mismatch_temp(1,:)-gamma*match_mismatch_temp(2,:));
map11_LA=map11_cell{ind_t};


        
        



