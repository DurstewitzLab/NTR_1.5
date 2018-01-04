function shuffled_data=shuffling(ordered_data,n_corr)
%Private, simple utility for cleaning the code. Shuffles rows 
%(and columns by uncommenting the corresponding code lines)
%of the ordered data. n_corr simply reffers to the number of consecutive
%time bins that has to be preserved during the shuffling.
[n_patterns,n_units]=size(ordered_data);
%Calculate groups, lower bound
n_groups=floor(n_patterns/n_corr);
%rand('state',0);

a=randperm(n_groups);%a permutation of the groups of 'n_corr' elements
shuffled_data=zeros(n_patterns,n_units);
for i=1:n_groups,
    b=randperm(n_units);
    shuffled_data(  (i-1)*n_corr+1 : i*n_corr , : )=...
        ordered_data( (a(:,i)-1)*n_corr+1 : a(:,i)*n_corr ,  b );
    %               shuffled_data(  (i-1)*n_corr+1 : i*n_corr , : )=...
    %                   ordered_data( (a(:,i)-1)*n_corr+1 : a(:,i)*n_corr ,:);
end
%Latest elements in 'shuffled_data_1' (if any) take the value of the last event stored in previous
%loop, i.e. shuffled_data_1(n_groups*n_corr,:
remaining=n_patterns-n_groups*n_corr;
if remaining>0
    %remaining=ceil(remaining);% Upper integer rounded
    shuffled_data(end-remaining+1:end,:)= ...
        repmat(shuffled_data(n_groups*n_corr,:),remaining,1);
end