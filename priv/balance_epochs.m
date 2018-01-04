function boot_lab_indexes=balance_epochs(epoch_labels,factor_disprop)
%Private, auxiliar file. Size of the longest epoch cannot exceed x "factor_disprop" the
%size of the shortest epoch.
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
%@ebalaguer: corrected april 012
%
%__________________________________________________________________________
%
warning('on','equilibrating_epochs:Eq')
if nargin<2,
    factor_disprop=10;
elseif factor_disprop<=0
    factor_disprop=inf;
    disp('Task-epochs retain their original size during analyses'),
end
m_lab=max(epoch_labels);min_lab=min(epoch_labels(epoch_labels>-1));l_lab=length(epoch_labels);
set_labs=[];n_bins_lab=[];
counter_labs=0;
for i=min_lab:m_lab
    if ~(all(i.*ones(l_lab,1)-epoch_labels))
        set_labs=[set_labs,i];
        counter_labs=counter_labs+1;
        n_bins_lab=[n_bins_lab,length(epoch_labels(epoch_labels==set_labs(counter_labs)))];
    end
end
n_lab=length(set_labs);
%
min_length=min(n_bins_lab); boot_lab_indexes=[];time_order=(1:l_lab)';
allowed_l=factor_disprop*min_length;
%
%Remove entire trajectories when needed
for i=1:n_lab,
    lab_indexes=time_order(epoch_labels==set_labs(i));%Real order corresponding to this label
    frac=n_bins_lab(i)/min_length;
    if (abs(frac)>factor_disprop),
        warning('equilibrating_epochs:Eq',['Epoch ',num2str(set_labs(i)),' length limited to x ',num2str(factor_disprop),' the smallest epoch']),
        %Removing the ones who exceed the criterion.
        %for that means, only entire trajectories can be randomly taken
        %out. Number of removed vectors have to be the same for each trial.
        removed=round( (n_bins_lab(i)-allowed_l) );
        lab_index_X_trial=time_order(epoch_labels==set_labs(i));
        %shuff_ind=shuff_data(lab_index_X_trial);%Randomly re-ordering but preserving consecutive vectors
        %shuff_ind=shuff_ind(removed:end);
        shuff_ind=lab_index_X_trial(1:(end-removed));
        remain_ind=sort(shuff_ind);
        boot_lab_indexes=[boot_lab_indexes;remain_ind];
    else
        boot_lab_indexes=[boot_lab_indexes;lab_indexes];
    end
end
boot_lab_indexes=sort(boot_lab_indexes);
