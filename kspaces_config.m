function conf=kspaces_config()
%v 4.0
%Configuration file.
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
%__________________________________________________________________________
%
%PLEASE SPECIFY IN THIS FILE THE CONFIGURATION PARAMETERS (1x1).
%
%Example of configuration:
%   (Multi)nomial kernel with expansion order 3.
%   Regularization 30% of the mean kernel matrix value.
%   Show multivariate discriminant statistics.
%   Displays activated of trial 1 and flow of trials 1-5 will be shown in the
%       maximum discriminating subspace
%   Full stat. analyses not reported
%For information about parameter settings please see comments below or 'readMeFirst.pdf'
%
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%1-Basic parameters, they have to be specified by the user
%
order=5;% (1x1) Natural number>0.  Degree of the product among neural responses across units.
%
is_shuffled_events=0;%(1x1) Natural. If >0, shuffles blocks as specified by the labels in �end-2" data column (task-epoch labels).
%
is_shuffled_within_events=0;%(1x1). Natural. If>0, preserves blocks as specified by the labels in �end-2" data column (task-epoch labels), but shuffles all time-bins within. 
%
trial_disp=5;%(1x1) Natural number indicating the trial to be displayed.
%                   Only one trial at a time is allowed to be displayed. 
%                   on a "3-D" trajectory plot. If <1, there will be no display at all.
%
make_video=5;%(1x1) Natural number indicating the trial to be displayed. 
%             If>0 a video of a single trajectory will be displayed
%             instead the static plot of the trial, regardless the value of
%             "trial_disp".
%
trials_flow_disp=[4];%(1x2) Range of trials "flows" to be displayed in three
%                   maximum discriminating (canonical) axis. If empty, no
%                   display will be produced.
%
%
%
%__________________________________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________
%
%
%
%
%2-Advanced parameters. No need to be modified unless default options are not enough.
%
regularization=.25;% (1x1) Dimensionality penalization (0,1), represtenting fraction of the average 
%                   value of the kernel matrix( see kfdmultiv.m for technical details)
%                   Suggested values: Lower values that default (e.g.
%                   0.148; 0.185) provide optimal cross-validations on
%                   weakly nonstationary time series(for instance see
%                   Balaguer-Ballester & Lapish et al. 2011)
%                   In contrast, larger values than default (e.g. 0.423; 0.396) typically provide best
%                   regularization results on non-stationary time series.

is_multiple_discriminant=1;% (1x1). Natural. If>0, prediction statistics based on a multiple discriminat
%             will be also reported (besides the average of all pairwise combinations of groups).
%
is_pairwise=0;%Same for the 2-class disc. If is_multiple_discriminant<0, the two-class will
%                   be performed anyway
%
xval_type=-1;% (1x1). Natural. %Type of Cross-validation that will be performed. 5  kinds of cross-validation are supported:
%                               0: Leave-one-out (standard method): each ith-trial is removed in turn then optimum discriminant
%                                                                directions are computed using all-but-the ith trial;
%                                                                this trial will be the ith-validation set
%                               1: Causal n-fold cross-validation 
%                               (respects causality): last j-trials are removed in each jth validation
%                                                     block. The remaining trials form the reference set, 
%                                                     which thus is smaller for increasing jth-validation blocks.
%                                                     There will be n=(m/2)-1 validation blocks (m=number of trials)
%                                2: Single-trial cross validation: each individual trial is a training sample and validation one
%                                                      in turn. There will be n=2*m validations blocks. This could be  used when 
%                                                      trials are not comparable for example when they do belong to different animals
%                               3: Single trial+causal. 
%                               
%                               <0: Only classification, trial by trial. Warning: use only once the regul.
%                                   constant have been optimized previously by cross-validation
%
full_stats_disp=1;%(1x1) Warning: Matlab Statistical Toolbox is needed for activating this option.
%                   If==1, command window information, some statistics as
%                   certain warnings will be displayed
%                   If > 1, all warnings will be displayed as well as multiple parametric and non-parametric tests 
%                   If>2, a full display of the procces will be shown
%
balance_trajec_length=30;%(1x1) This parameters controls the disbalanced in epochs size:
%                 If>0, the number of time bins of the longest
%                 task-epochs is, at maximum "x balance_task_epochs" of the shortest
%                 one. Useful when classes have very different sizes.
%                 In addition, if this parameter is positive, chance priors will 
%                 be automatically used in the discriminant analyses (for
%                 avoiding size-biases between classes).            
%  
do_dcm=1;%(1x1) If>0, performs a delay-coordinate map expansion before kernel expansion.
%           Not compelling in multivariate recordings unless all
%           units responses are often the same for two different
%           task-epochs. In that situation, trajectories would better be
%           disambiguated by a lag-expansion.
%           Note: Regardless this parameter value,
%           a delay-coordinate map will be used in the trajectory visualization depending 
%           only on the value of the next parameter.           
%
fixed_lag=30;%(1x1) Note: Parameter only used during 3D trajectory visualization pruposes, 
%           not for the high-dim statistical analyses neither for the flow display. 
%           Natural number, indicating the maximum lag (in time
%           bins units) which should have every axes in a delay-coordinate map. If <=0, lags will be
%           obtained as the average/maximum lag across the "first" minimum between cross-correlation 
%           functions (see code below in dcm.m for details).
%   
incl_subop_lags=1;%(1x1) Note: Parameter only used for the high-dim statistical analyses & flow display,
%                 not during 3D trajectory display. 
%                 Binary, indicating if all lags below the
%                 optimum/maximum one will be included as a new axes or not. 
%
balance_task_epochs=30;%(1x1) This parameters controls the disbalanced in epochs size:
%                 If>0, the number of time bins of the longest
%                 task-epochs is, at maximum "x balance_task_epochs" of the shortest
%                 one.
%
lost_frac=.1;%(1x1) Note: Parameter only used for the high-dim statistical analyses & flow display,
%             not during 3D trajectory display. 
%             Maximum percentage of lags that can be lost during a DCM embedding.
%
do_record=0; %(1x1) Record a video of the trajectory display in mp4 format
%
video_quality=100;% (1x1) Quality of the recorded video (if any) in (0-100%)
%
is_gauss_kernel=0;%(1x1) WARNING: Not implemented, found not necessary. Value will be ignored & can be ommited]
%             If >0, infinite dimensional gaussian kernel will be used (not interpretable as
%             spike matrix but translation-invariant. The value of this
%             configuration, when >=0, specifies standard deviation (in units of bins)
%             used in the kernel. if is_expo>=0, polynomial kernel will be not used
%             then, thus "order" parameter will be irrelevant.
%__________________________________________________________________________

conf=setup_config(order,regularization,...
    is_shuffled_events,is_shuffled_within_events,...
trial_disp,make_video,trials_flow_disp,is_multiple_discriminant,is_pairwise,xval_type,...    
full_stats_disp,balance_trajec_length,do_dcm,fixed_lag,incl_subop_lags,balance_task_epochs,...
lost_frac,do_record,video_quality,is_gauss_kernel);
end


