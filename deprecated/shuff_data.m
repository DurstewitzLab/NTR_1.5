function shuff=shuff_data(data,shuff_epoch,shuff_within_epoch)
%Shuffles across epochs, either preserving
%each individual task-epoch temporal contingency
%
%Based on: Balaguer-Ballester E, Lapish C, Seamans JK and Durstewitz D 2011 
% "Attracting Dynamics of Frontal Cortex Populations during Memory Guided Decision-Making". 
% PLoS Comput. Biol.
%by BCCN Heidelberg-Mannheim.
%
%Updates will be provided in this address:
%http://www.bccn-heidelberg-mannheim.de/
%emili.balaguer@zi-mannheim.de
%Contact adresses: emili.balaguer@zi-mannheim.de
%                  daniel.durstewitz@zi-mannheim.de
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%
%Inputs:
%   data (number of time bins x (dimensionality of neural responses+3) ), and with the following
%                                                                         structure:
%   (note: Dimensionality of neural responses refers e.g. to the number of simultaneously
%    recorded "neurons" and/or delayed version of them)
%           -Columns "1:end-3" of "data" matrix:  Must contain neural responses values over
%                                       time.
%           -Column "end-1" of "data" matrix: Must contain natural numbers, they are
%                                      alternative labelling used only for
%                                      display (number of
%                                      labels smaller than 8). If all "epochs"
%                                      are to be displayed, this column has
%                                      to be a copy of "end-2" one.
%                                      NOTE: these labels will determine the groups
%                                      displayed in the smoothed
%                                      trajectory (upper plots). 
%           -Column "end-2" of "data" matrix: Must contain natural numbers, labelling
%                                      the different stimulus or behavioral "epochs" in which
%                                      the experimentalist segments the
%                                      task. "-1" encodes "no-labelled"
%                                      time-bins. NOTE: used for colouring
%                                      arrows in lower plots.
%           -Last column of "data" matrix: Must contain natural numbers, labelling the
%                                      different trials of the task.
%
%   shuff_epoch (1x1) Real. If >0, shuffles blocks as specified by the first column in labels.
%   shuff_within_epoch (1x1) Real. If>0, preserves blocks as specified by the first column
%                               within labels and shuffles withi those
%                               blocks.
%
%Outputs:
%shuffled (number of time bins x (dimensionality of neural responses+3)
%       and with exactly the same structure as the input "data" matrix.
%__________________________________________________________________________

%I. Preliminaries
stream = RandStream('mt19937ar');
RandStream.setDefaultStream(stream);
%
if nargin<3, shuff_within_epoch=0; end
if nargin<2, shuff_epoch=1; end
trajec=[];%A single trajectory of a given task-epoch
trajec_acum={};%Total number of trajectories is unknown. Each cell is a trajectory.
trajec_durations=[];%In time bins
trajec_length=0;%Counter for the number of different trajectories, of any kind
max_trajec_length=30;%Maximum allowed trajectory length, in bins (hard-coded)
%
%Recovering labels of trajectories and activity vectors data. Remember that last column
%contains the  number of trial, previous-to-last class-labels only for trajectory plotting,
%whereas the second-to-last contains the task-epoch labels which will be analysed in here.
%
%However this funcion can be used for block-shuffling single vector of time indexes. 
[~,columns]=size(data);
if columns>3
    index_trajec=data(:,end-2);
    activ=data(:,1:end-3);
elseif columns==2
    index_trajec=data(:,2);
    activ=data(:,1);
elseif columns==1
    index_trajec=data;
    activ=index_trajec;
else
    error('Input data can be a vector of tiem indexes or a data matrix of neural responses with either 1 or 3 final labels columns')
end   

n_bins=length(activ(:,1));
%
%II. Create a cell array, each cell is a trajectory
for i=1:n_bins,
    if ( (i==1) || ( (index_trajec(i)==index_trajec(i-1)) ) ),%Those two vectors belong to the
        %                                                      same trajectory.
        %
        if (i==n_bins)||(trajec_length==max_trajec_length),
            
            if (shuff_within_epoch>0)
                warning('shuff_data:withinTrShuf','Data altered by randomly shuffling within individual task-epoch trajectories'),
                trajec=shuffling(trajec,1);%Activate to shuffle also within trajectories.
                warning('off','shuff_data:withinTrShuf')
            end
            trajec_durations=[trajec_durations,length(trajec(:,1))];%Stores epoch durations for step 'II'.
            trajec_acum=[trajec_acum,trajec];
            trajec_length=0; trajec=activ(i,:);
        else
            trajec=[trajec;activ(i,:)]; trajec_length=trajec_length+1;
        end
    else
        if (shuff_within_epoch>0)
            warning('shuff_data:withinTrShuf','*Data altered by randomly shuffling within individual task-epoch trajectories*'), %
            trajec=shuffling(trajec,1);%Activate to shuffle also within trajectories.
            warning('off','all')
        end
        trajec_durations=[trajec_durations,length(trajec(:,1))];%Stores epoch durations for step 'II'.
        trajec_acum=[trajec_acum,trajec];
        trajec=activ(i,:);
    end
end
%adding the last vector to the last trajectory
trajec_durations(end)=trajec_durations(end)+1;
%Add the last one
trajec_acum=[trajec_acum(1:end-1),[trajec_acum{end};activ(n_bins,:)]];

%
%III. Shuffling by preserving blocks of trajectories
n_trajec=length(trajec_acum);
shuffled_activ=[];%Matrix, of exactly the same dimension as "activ" and containing the shuffled
%                 trajectories.
permuted=randperm(n_trajec);
for j=1:n_trajec,%Loop over all trajectories
    if shuff_epoch>0
        warning('shuff_data:acrossTrShuf','*Data altered by randomly shuffling entire blocks of task-epoch trajectories*'),
        warning('off','shuff_data:acrossTrShuf')
        shuffled_trajec=trajec_acum{permuted(j)};
    else
        shuffled_trajec=trajec_acum{j};
    end
    shuffled_activ=[shuffled_activ;shuffled_trajec];
end
%
%Adding the three last columns of the original (input) matrix, containing labels
if columns>3
    shuff=[shuffled_activ,data(:,end-2:end)];
elseif columns==2
    shuff=[shuffled_activ,data(:,end)];
elseif columns==1
    shuff=shuffled_activ;
end
warning('on','all')

