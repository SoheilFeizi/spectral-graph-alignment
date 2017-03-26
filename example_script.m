% An example script to run EigenAlign and LowRankAlign
% please see individual functions for more details

clc
close all
clear all

% generating example graphs

n=50;
rng(1);
N1=rand(n)>0.8;
N1=triu(N1);
N1=N1+N1';
N1=N1~=0;

for i=1:n
    N1(i,i)=0;
end

map_true=zeros(n);
for i=1:n
    map_true(i,n-i+1)=1;
end

N2=N1;
% making reverse ordering
N2=N2(n:-1:1,n:-1:1);

%****************
gamma=0;
%running EigenAlign using default parameters
disp('running EigenAlign')
map=ones(n);
map11_EA=EigenAlign(N1,N2,map,gamma);

figure
spy(map_true)
title('true mapping')

figure
spy(map11_EA{1})
title('inferred mapping-EigenAlign')

%************************
%running lowrank align
k_rank=2;% number of top eigenvectors considered

map11_LA=LowRank_Align(N1,N2,map,k_rank,gamma);

figure
spy(map11_LA)
title('inferred mapping-LowRank Align')

