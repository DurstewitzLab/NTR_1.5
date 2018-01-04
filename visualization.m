function visualization(data,order,trial_disp,fixed_lag,do_dcm,make_video,trials_flow_disp,...
    balance_task_epochs,regul,incl_subop_lags,lost_frac,do_record,videoQuality)
%v 2
%Displays a "single-trial" (short time series) using kernel-PCA (Sch�lkopf et al., 1998)
%and "flow" or "velocity vectors" of multiple-trials using kernel-FDA (Mika et al., 2000).
%
%Orders #1 and higher "order" will be displayed in a 2x2 Figure.
%It is only an orientative display, not precisely reflecting high-dimensional 
%detailed statistics.
%
%This code is not commercial and thus is not guaranteed but is used on regular basis in 
%multiple datasets.
%Based on: Balaguer-Ballester E, Lapish C, Seamans JK and Durstewitz D 2011
% "Attracting Dynamics of Frontal Cortex Populations during Memory Guided Decision-Making".
% PLoS Comput. Biol.
%Distributed accourding to the General Public License GNU 3.0
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
%Inputs:
%   data (number of time bins x (dimensionality of neural responses+3) ), and with the following
%                                                                         structure:
%
%   (note: dimensionality of neural responses refers e.g. to the number of simultaneously
%    recorded "neurons" and/or delayed version of th3m)
%
%           -Columns "1:end-3" of "data" matrix:  Must contain activity values over
%                                       time.
%           -Column "end-1" of "data" matrix: Must contain natural numbers>0, they are
%                                      alternative labelling used only for
%                                      display (number of
%                                      labels smaller than 8). If all "epochs"
%                                      are to be displayed, this column has
%                                      to be a copy of "end-2" one.
%                                      NOTE: these labels will determine the groups
%                                      displayed in the smoothed
%                                      trajectory (upper plots). 
%           -Column "end-2" of "data" matrix: Must contain natural numbers>0, labelling
%                                      the different stimulus or behavioral "epochs" in which
%                                      the experimentalist segments the
%                                      task. "-1" encodes "no-labelled"
%                                      time-bins. NOTE: used for colouring
%                                      arrows in lower plots.
%           -Last column of "data" matrix: Must contain natural numbers>0, labelling the
%                                      different trials of the task.
%   order (1x1) Natural number>0. Degree of the product among spike densities across units.
%
%   trial_disp (1x1) Natural."3D" smoothed trajectories display.
%                      Only one trial at a time is allowed to be displayed.
%   fixed_lag (1x1) Natural number, indicating the maximum lag (in time
%           bins units) which should have every axes. If <=0, lags will be
%           obtained as the average/maximum lag across the "first" minimum between cross-correlation 
%           functions (see code below in dcm.m for details).
%
%   do_dcm (1x1) Important note: Only relevant for the flow display, not to the
%       trajectory display (in which there is always a dcm performed). Performs a delay-coordinate map expansion before kernel
%       expansion. 
%
%   make_video 1x1) Natural number indicating the trial to be displayed. 
%             If>0 a video of a single trajectory will be displayed
%             instead the static plot of the trial, regardless the value of
%             "trial_disp".
%
%   balance_task_epochs (1x1) This parameters controls the disbalanced in epochs size:
%                 If>0, the number of time bins of the longest
%                 task-epochs is, at maximum "x balance_task_epochs" of the shortest
%                 one.
%
%   trials_flow_disp (1x2) First and last trials "flow" to be displayed in three
%                   maximum discriminating (canonical) axis. If<1, no
%                   display will be produced.

%   regul (1x1) Regularization constant, only fo flow display
%__________________________________________________________________________
figure
set(gcf,'Color',[1,1,1],'name',...
     ['Single-trial complete trajectory and vector field for the activity space and expanded space (order=',num2str(order),',max. lag ',num2str(fixed_lag),'bins)']),
%
activ=data(:,1:end-3);
trial_labels=data(:,end);
phase_labels=data(:,end-1);%remaining: "phase_labels" are those used onyl for trajectory display.
epoch_labels=data(:,end-2);
%
%I. SMOOTHED TRAJECTORIES OF A SINGLE TRIAL
disp('_____________________________')
disp('*PERFORMING VISUALIZATIONS...*')
if (trial_disp>0)
    %
    %I.1. Checking that such a trial exists
    if trial_disp>max(trial_labels)
        warning('visualization:InvalidTr','This trial does not exists, no trajectory plot will be produced')
        subplot(221), text(0.5,0.5,'This trial does not exists')
        subplot(222), text(0.5,0.5,'This trial does not exists')
    else
    %I.2 Select the activity from this trial
        activ_tr=activ((trial_labels==trial_disp),:);
        phase_labels_tr=phase_labels((trial_labels==trial_disp),:);
        %
        %I.3 Reduced data for activity product order=1
        warning('off','all')
        del_activ=dcm(activ_tr,fixed_lag);
        indexes=balance_epochs(phase_labels_tr(1:length(del_activ(:,1))),balance_task_epochs);%Returns the indexes of the radomly selected labels
        del_activ=del_activ(indexes,:);
        reduced_data=kpca(activ_tr);
        warning('on','all')
        display_trial(reduced_data,phase_labels_tr(1:length(reduced_data(:,1))),trial_disp,1,1);
        %
        %I.4 Using a delay-coordinate map and data reduction for the selected order.
        warning('off','all')
        reduced_data=kpca(del_activ,order);
        warning('on','all')
        %
        %Display a video or a plot for the selected activity expansion
        %order.
        if make_video>0
            video_trial(reduced_data,phase_labels_tr(1:length(reduced_data(:,1))),trial_disp,order,2);
            %
            if do_record>0
                %If recording a trial, create a new fig. and close it
                figure
                video_trial(reduced_data,phase_labels_tr(1:length(reduced_data(:,1))),trial_disp,order,-1,do_record,videoQuality);
                close (gcf);
            end
        else
            display_trial(reduced_data,phase_labels_tr(1:length(reduced_data(:,1))),trial_disp,order,2);
        end
    end
    disp('...3D trajectory displayed')
end

%II. TRAJECTORY FLOW

if length(trials_flow_disp)==2,
    %
    %II.1. Checking that there are not too many trials to display
    if (max(trials_flow_disp)>max(trial_labels))||(min(trials_flow_disp<1)),
        warning('visualization:InvalidTrRange','Invalid trial range selected, no "flow" plot will be produced')
        subplot(223), text(0.5,0.5,'Invalid trials range')
        subplot(224), text(0.5,0.5,'Invalid trials range')
    else
        %
        %II.2 Select the activity from several epochs
        activ_tr=activ(  ((trial_labels>=trials_flow_disp(1))&(trial_labels<=trials_flow_disp(2))),:);
        labels=epoch_labels(  ((trial_labels>=trials_flow_disp(1))&(trial_labels<=trials_flow_disp(2)))  );
        %
        %Use only this map if requested for the statistics
        del_activ=activ_tr;
        if do_dcm>0
            warning('off','all'),
            %Use or not fixed lag may improve visualization
            do_display=1;
            del_activ=dcm(activ_tr,0,incl_subop_lags,do_display,lost_frac);
            %Note: las argument indicates "smoothing/non-smootihg" of
            %trajetories.
            warning('on','all'),
        end
        l=length(del_activ(:,1));
        labels=labels(1:l,:);
        %
        %II.3 Reduced data in the maximum discriminating subspace for
        %activity product order=1 and for the selected "order" (using multi-class regularized
        %Fisher Discriminant in high dimensions)
        %
        warning('off','all'),
        indexes=balance_epochs(labels,balance_task_epochs);%Returns the indexes of the radomly selected labels
        del_activ=del_activ(indexes,:);labels=labels(indexes,:);
        reduced_data=kfd_multiv(activ_tr(indexes,:),labels,1,regul);
        warning('on','all'),
        %display_flow(reduced_data,labels,trials_flow_disp,1,3);%deprecated
        display_flow2(reduced_data,labels,trials_flow_disp,1,3);
        warning('off','all'),
        reduced_data=kfd_multiv(del_activ,labels,order,regul);
        warning('on','all'),
        %display_flow(reduced_data,labels,trials_flow_disp,order,4);%deprecated
        display_flow2(reduced_data,labels,trials_flow_disp,order,4);
        pause(5),
    end
end
%Example of figure name, should you want to use it
set(gcf,'Color',[1,1,1],'name',...
     ['Single-trial complete trajectory (order=',num2str(order),',max. lag=',...
     num2str(fixed_lag),'bins, regul ',num2str(regul),')']),
