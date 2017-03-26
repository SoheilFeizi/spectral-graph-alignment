function [num_matches,num_mismatches]=eval_rate_match_mismatch2(N1,N2,XX)


N22=XX'*N1*XX;
N11=XX*N2*XX';

num_matches=sum(sum(N22.*N2))/2;
num_mismatches=sum(sum(N22.*(1-N2)))/2+sum(sum(N11.*(1-N1)))/2;