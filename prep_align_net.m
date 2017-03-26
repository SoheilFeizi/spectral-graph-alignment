function A=prep_align_net(N,M,P)

% constructing alignment network
% N and M are input binary networks (square matrices), possibly asymetric
% N and M do not have self loops
% P is a mapping across networks

n=size(N,1);
m=size(M,1);

%**************
% binarize networks if not binary
N=(N>0);
M=(M>0);
P=(P>0);

%**************
% networks do not have self loops
% N=N.*(1-eye(n));
% M=M.*(1-eye(m));

for i=1:n
    N(i,i)=0;
end

for i=1:m
    M(i,i)=0;
end

%**************
[x,y]=find(P~=0);
z=[x,y]; % list of matches
nz=size(z,1);

%**************
% use a sparse version if networks are very large
% or use LowRank Align if too slow
%nz
if nz>50000
    disp('alignment graph is too large, switching to sparse mode')
    
    disp('if too slow, use LowRank_Align.m')
    A=sparse(nz,nz,0,nz,nz);
else
    A=zeros(nz,nz);
end;
for k1=1:nz-1
    for k2=(k1+1):nz
        
        x1=x(k1);
        x2=x(k2);
        
        y1=y(k1);
        y2=y(k2);
        
        if min(N(x1,x2)*M(y1,y2)+N(x2,x1)*M(y2,y1),1)==1;
            % matches
            A(k1,k2)=1;
            A(k2,k1)=1;
            
        elseif N(x1,x2)~=0 || N(x2,x1)~=0 || M(y1,y2)~=0 || M(y2,y1)~=0
            % mis-matches
            A(k1,k2)=-1;
            A(k2,k1)=-1;
        end
    end
end

%*************
A=sparse(A);
