function d_data=lagged_axes(data,lags)
%A minimal embedding can be acheived just by substituting original axis  with
%the optimum delayed ones ("whitening by using delays"). It migth be useful
%in cases where axes are correlated but none can be really discardedand
%dimensionality do not have to be increased
%
%This code is not commercial and thus is not guaranteed but is used on regular basis in 
%multiple datasets.
%Based on: Balaguer-Ballester E, Lapish C, Seamans JK and Durstewitz D 2011
% "Attracting Dynamics of Frontal Cortex Populations during Memory Guided Decision-Making".
% PLoS Comput. Biol.
%Distributed accourding to the General Public LIcense GNU 3.0
%http://www.gnu.org/copyleft/gpl.html
%
%by Emili Balaguer-Ballester.
%Contact adresses: eb-ballester@bournemouth.ac.uk
%                  daniel.durstewitz@zi-mannheim.de
%                  emili.balaguer@uv.es
%
%
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%
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

    