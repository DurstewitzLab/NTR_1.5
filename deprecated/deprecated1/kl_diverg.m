function KL=kl_diverg(observed,expected)
%march 011
%Naive information-theoretic test of differences between two pdf's
%Discretized KL(distribution1,distribution2) is defined as
%sum for all 'i' of:  (distribution1(i)*log(distribution1(i)/distribution2(i)).
%Then is zero if distribution1(i)=distribution2(i) for each 'i'
%Inputs: histograms of both distributions (multi-dim).
%Output: (1x1) scalar contaning the euclidean length of the KL divergence vector. 
warning('off','kl_diverg:zeroCounts'),
%Setting the same size of both histograms
[n1,columns_o]=size(observed);[n2,columns_e]=size(expected);
if (n1~=n2), warning(20,'Both histograms forced to have the same length'),
end
 if n1>n2, n=n2; else n=n1; end
if (columns_o<columns_e), expected=expected(:,1:columns_o); columns=columns_o;
else  observed=observed(:,1:columns_e); columns=columns_e; end
%
%Shouts if any count is zero
minimum_q=min(min(observed));minimum_p=min(min(expected));
if (minimum_q<=0)||(minimum_p<=0),
    warning('kl_diverg:zeroCounts','Some counts=0 replaced by the smallest count value, then histograms are normalized')
    min_o=min(observed(observed>0)); observed(observed==0)=min_o; 
    min_e=min(expected(expected>0)); expected(expected==0)=min_e; 
end
%
%(re-) normalization for producing proper pdfs per dimension
%i.e. separately nromalizing each column
observed_dist=observed./repmat(sum(observed,1),n,1);
expected_dist=expected./repmat(sum(expected,1),n,1);
KL=norm(observed_dist'*log(observed_dist./expected_dist));
