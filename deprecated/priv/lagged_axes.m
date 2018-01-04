function d_data=lagged_axes(data,lags)
%A minimal embedding can be acheived just by substituting original axis  with
%the optimum delayed ones ("whitening by using delays"). It migth be useful
%in cases where axes are correlated but none can be really discardedand
%dimensionality do not have to be increased
[n_patterns,n_channels]=size(data);
d_data=zeros(n_patterns,n_channels+1); %Includes the undelayed first axis and the 2-delayed one (in lags(1)).                                                                                     after that (lags(2)...lags(n_channels+1)
%                                      % After that (lags(2)...lags(n_channels+1)contains a 1-delay for each lag.
%
d_data(:,1)=data(1:n_patterns,1);%1st column of 'd_data' matrix is the first time series undelayed.
lost_vectors=zeros(n_channels,1);
for j=2:n_channels+1,
    lost_vectors(j-1)=lags(j-1)-1;
    d_data(1:n_patterns-lost_vectors(j-1),j)=data(   lags(j-1):n_patterns,j-1    );
end
d_data=d_data(1:n_patterns-max(lost_vectors),:);

    