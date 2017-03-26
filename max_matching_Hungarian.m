function XX=max_matching_Hungarian(X)

% X is the input network

X_r2=-X;
X_r2(X==0)=inf;
[XX,~] = Hungarian(X_r2);
