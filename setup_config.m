function conf=setup_config(order,regularization,...
    is_shuffled_events,is_shuffled_within_events,...
    trial_disp,make_video,trials_flow_disp,is_multiple_discriminant,is_pairwise,xval_type,...
    full_stats_disp,balance_trajec_length,do_dcm,fixed_lag,incl_subop_lags,balance_task_epochs,...
    lost_frac,do_record,video_quality,is_gauss_kernel)
%v 4.0
%Auxiliar configuration file. 
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
%For information  about parameter settings please see 'ReadMe.txt'
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%@ebballester: Error in attribute 9, it was the number 10
%
%__________________________________________________________________________
%
%Inputs:
%
%   I.BASIC CONFIGURATION parameters
%
%         order (1x1). Natural number>0. Degree of the product among spike densities across units.
%
%         regularization (1x1). Real. Penalizes the complexity of the Fisher Discriminant Criterion in units of
%             percent of the mean kernel matrix. If this value is
%             <0, regularization will be optimized such that the
%             leave-one-out error is minimized [WARNING: This automatisation is not implemented in
%             v4.0b. Will thrown an error]
%
%         is_shuffled_events (1x1). Real. If >0, shuffles blocks as specified by the first column in labels.
%
%         is_shuffled_within_events (1x1). Real. If>0, preserves blocks as specified by the first column in labels
%                 and shuffles within blocks.
%
%         trial_disp (1x1). Natural number. Only one trial at a time is allowed to be displayed.
%             on a "3-D" trajectory plot. If <1, there will be no display at all.
%
%         make_video (1x1). Natural, if>0 a video of a single trajectory will be
%              displayed (regardless the value of "trial_disp").
%
%         trials_flow_disp (1x2).  Number of trials "flows" to be displayed in three
%                   maximum discriminating (canonical) axis. If =[], no
%                   flow display will be produced.
%
%_  _  _  _  _  _ _  _  _  _  _  _ _  _  _  _  _  _ _  _  _  _  _  _ _  _ _
%_  _  _  _  _  _ _  _  _  _  _  _ _  _  _  _  _  _ _  _  _  _  _  _ _  _ _
%
%   II.ADVANCED CONFIGURATION parameters
%
%       is_multiple_discriminant (1x1). If>0, prediction statistics based on a multiple discriminant
%             will be also reported (besides the average of all pairwise combinations of groups).
%
%           xval_type (1x1). Type of Cross-validation that will be performed. Three  kinds of cross-validation are supported:
%                               0: Leave-one-out (standard method): each ith-trial is removed in turn then optimum discriminant
%                                                                directions are computed using all-but-the ith trial;
%                                                                this trial will be the ith-validation set
%                               1: Causal n-fold cross-validation
%                               (respects causality): last j-trials are removed in each jth validation
%                                                     block. The remaining trials form the reference set,
%                                                     which thus is smaller for increasing jth-validation blocks.
%                                                     There will be n=(m/2)-1 validation blocks (m=number of trials)
%                                >1: Single-trial cross validation: each individual trial is a training sample and validation one
%                                                      in turn. There will be n=2*m validations blocks. This could be  used when
%                                                      trials are not comparable for example when they do belong to different animals
%
%       full_stats_disp (1x1). Warning: Matlab Statistical Toolbox is needed for activating this option.
%           If>0, multiple parametric and non-parametric tests will be displayed in the command window.
%
%       balance_trajec_length (1x1) = If >0, the length of any class-specific-trajectory is
%           upper-bounded by the smallest of the mean trajectory length (across all groups) and
%           this number of bins. Default value is 0.
%           Useful when classes have very different sizes.
%
%       do_dcm (1x1). Performs a delay-coordinate map expansion before kernel expansion.
%           Not used by default in multivariate recordings unless all
%           responses are often the same for two different
%           taks-epochs. In that situation, trajectories should be
%           disambiguated by a lag-expansion.
%           In addittion, of this parameter is positive, chance priors will
%           be automatically used in the discriminant analyses (for
%           avoiding size-biases between classes).
%
%       fixed_lag (1x1). Natural number, indicating the maximum lag (in time
%           bins units) which should have every axes. This is useful just for
%           smoothing the trajectories on 3D displays. If <=0, lags will be
%           obtained as the average between the "first" minimum between cross-correlation
%           functions (see code below for more details).
%           It is also used for smooth 3D trajectory visualization, even if
%           do_dcm=0.
%
%       incl_subop_lags (1x1). Binary, indicating if all lags below the
%           optimum/maximum one will be included as a
%           new axes or not. It is also used for smooth visualization, even if
%           do_dcm=0.
%
%       balance_task_epochs.(1x1). This parameter controls the disbalanced in epochs size:
%                      If>0, the number of time bins of the longest
%                      task-epochs is, at maximum "x balance_task_epochs" of the shortest
%                       one.
%
%       is_gauss_kernel (1x1) [WARNING: Not implemented in v4.0b. Will be ignored]
%             If >0, infinite dimensional gaussian kernel will be used (not interpretable as
%             spike matrix but translation-invariant. The value of this
%             configuration, when >=0, specifies standard deviation (in units of bins)
%             used in the kernel. if is_expo>=0, polynomial kernel will be not used
%             then, thus "order" parameter will be irrelevant.
%Outputs:
%       conf (1x13) Structure in which each field contains upper parameters.
%__________________________________________________________________________

disp('_________________________________')
disp('*CURRENT CONFIGURATION*')
%__________________________________________________________________________
%I-BASIC CONFIGURATION
%
if (nargin<1)||(order<1)||(order>10)
    warning('setup_config:highOrder','Too high/low multinomial order. Set to "1"')
    order=1;%Natural number>0.  Degree of the product among spike densities across units.
else
    disp(['   Multinomial expansion order ',num2str(order)]),
end
if (nargin<2)||(regularization<0) %||(regularization>10)
    %Positive and real. Penalizes the complexity of the Fisher
    %Discriminant Criterion.
    warning('setup_config:noReg','Unfeasible regularization parameter')
    error('Regularization automatically optimized is not yet implemented. Please choose a positive parameter'),
else
    disp(['   Regularization is ',num2str(regularization*100),'% of the mean kernel matrix']),
end
if (nargin<3)||(is_shuffled_events<=0)
    is_shuffled_events=0;
else
    disp('   Shuffle data within task-epochs');%Suffles blocks as specified by the first column in labels
end
if (nargin<4)||(is_shuffled_within_events<=0)
    is_shuffled_within_events=0;
else
    disp('   Shuffle data within task-epochs');%Preserves blocks as specified by the first column in labels and shuffles within them
end
if (nargin<5)||((trial_disp<=0))
    %Different examples of displays
    if (make_video<=0)
        warning('setup_config:NoTrDisp','No trial trajectory displayed'),
    end
else
    disp(['   Display trial ',num2str(trial_disp),' smoothed trajectory']);
end
if (nargin<6)||( make_video<=0);
    %No Movie.
else
    %The trial indicated in here will be the true one to be
    %displayed, regardless the value of "trial_disp".
    %%No simultaneous video-plot display.
    trial_disp=make_video;
    disp(['   Make video of trial ',num2str(trial_disp),' trajectory']);
end
if (nargin<7)||(length(trials_flow_disp)<1)||(trials_flow_disp(1)<=0)
    warning('setup_config:noFlowDsip','No task-epochs flow displayed');
else
    if length(trials_flow_disp)<2,trials_flow_disp(2)=trials_flow_disp(1); end
    if (trials_flow_disp(2)<=0)
        warning('setup_config:noValid','Upper trial limit not valid, only one trial displayed');
        trials_flow_disp(2)=trials_flow_disp(1);
    end
    disp(['   Displaying flow of trials ',num2str(trials_flow_disp(1)),...
        ' to ',num2str(trials_flow_disp(2))]);
end
%
%__________________________________________________________________________
%II-ADVANCED CONFIGURATION
%
if (nargin<8)||(is_multiple_discriminant<=0)
    is_multiple_discriminant=0;
    disp('   No Multiple discriminant-based statistics will be reported');
else
    disp('   Multiple discriminant-based statistics will be reported');
    %Compute the prediction error based on a multiple discriminat
end
if ( (nargin<9)||(is_pairwise<=0) ),
    is_pairwise=0;
     disp('   No two-class discriminant-based statistics will be reported');
else
    disp('   Two-class discriminant-based statistics will be reported');
    %Compute the prediction error based on a two-class discriminat
    %besides statistic across all pairwise combinations of groups.
end
disp('   Cross-Validation variant selected:'),
if (nargin<10) ||(xval_type==0)
    xval_type=0;
    warning('setup_config:noncausal','Non-causal predictions (standant Leave-one-out cross-validation'),
elseif (xval_type==1),
    disp('      Causal predictions reported (please see more info. in readme.txt)');
    xval_type=1;
elseif (xval_type==2),
    disp('      Non-causal, leave-one-out prediction in which both reference and validation sets are made of a single-trial (please see more info. in readme.txt)');
    xval_type=2;
elseif (xval_type<0),
    disp('      Only penalized classification '),
    warning('setup_config:no_cross','No cross validation performed: please make sure that regularization is correct'),
    xval_type=-1;
else
    disp('      Causal, leave-one-out prediction in which both reference and validation sets are made of a single-trial (please see more info. in readme.txt)');
    xval_type=3;
end
if (nargin<11)||(full_stats_disp<=0)
    %Do not perform any delay-coordinate map by default.
    disp('   Ommiting command window information and full statistical analyses')
else
    warning('setup_config:full_stats','Adittional command window information and full statistical analyses may be shown')
    disp('   Please set full_stats_disp=0 for simpler command window reports')
end
if (nargin<12)
    %Do balance trajectories by default
    balance_trajec_length=30;
    disp(['   Task-epoch specific trajectories length constrained to ',num2str(balance_trajec_length),...
        ' time bins']),
    disp('   Equal prior probabilities used in the analyses'),
elseif (balance_trajec_length<=0)
    warning('setup_config:TrUnc','Task-epoch specific trajectories length unconstrained and real prior probabilities used for clasification')
    disp('   Analyses could be biased by strongly different group sizes'),
else
    disp(['   Task-epoch specific trajectories length constrained to ',num2str(balance_trajec_length),...
        ' time bins']),
    disp('   Equal prior probabilities used in the analyses'),
    %     warning(105,['Task-epoch specific trajectories length constrained to ',num2str(balance_trajec_length),...
    %         'standart deviations over the smallest mean length for all taks-epochs'])
end
if (nargin<13)||(do_dcm<=0)
    %Do not perform any delay-coordinate map by default.
    disp('   Delay coordinate map not used for statistical analyses (only for displays)'),
else
    warning('setup_config:Dcm','Delay-coordinate map used for statistical analyses'),
end
if (nargin<14)
    fixed_lag=30;
    disp(['   Only for visualization the maximum lag is set to ',num2str(fixed_lag),' time bins']),
elseif (fixed_lag<=0)&&(do_dcm>0)
    warning('setup_config:OptimLags','Optimum lags optimized via cross-correlation between dimensions')
elseif (fixed_lag>0)
    disp(['   Only for visualization the maximum lag set to ',num2str(fixed_lag),' time bins']),
end
if (nargin<15)
    incl_subop_lags=1;%Include suboptim. lags by default
    if (do_dcm>0), warning('setup_config:lowLags','Lags below the maximum included in the statistical analyses'), end
elseif (incl_subop_lags<=0)
    if (do_dcm>0), disp('   Sub-optimal lags not included in the statistical analyses'), end
end
if (nargin<16)||(balance_task_epochs<=0);
    balance_task_epochs=0;%Task-epochs are not balanced
    disp('   Task-epochs orignal sizes preserved'),
else
    disp(['   Task-epochs maximum size cannot exceed x ',num2str(balance_task_epochs),' times the length of the shortest one']),
end
if (nargin<17)||(lost_frac<=0);
    lost_frac=0.003;
else
        disp(['   Maximum ',num2str(lost_frac*100),'% of lags lost during embedding for discriminant analyses & flow visualizazion']),
end
if (nargin<18)||(do_record<=0);
    do_record=0;
else
     disp('   .avi video recorded and stored in the current directory'),
end
if (nargin<19)||(video_quality<=0);
    video_quality=50;
else
     disp(['   Movie quality is ',num2str(video_quality),'%. The largest, the bigger is the video file stored']),
end
if (nargin<20)||(is_gauss_kernel<=0);
    is_gauss_kernel=0;%Different examples of displays
else
    %     disp(['Warning: Exponential kernel used instead of polynomial. ',...
    %         ' Thus "order" is irrelevant), witdth ',num2str(is_exponential_kernel)]);
    error('Gaussian kernel not implemented. Only multinomial available')
end
%
%Storing in a config object
conf.order=order;
conf.regul=regularization;
conf.shuff_epoch=is_shuffled_events;
conf.shuff_within_epoch=is_shuffled_within_events;
conf.trial_disp=trial_disp;
conf.make_video=make_video;
conf.trials_flow_disp=trials_flow_disp;
conf.multi=is_multiple_discriminant;
conf.pairwise=is_pairwise;
conf.xval_type=xval_type;
conf.full_stats_disp=full_stats_disp;
conf.balance_trajec_length=balance_trajec_length;
conf.do_dcm=do_dcm;
conf.fixed_lag=fixed_lag;
conf.incl_subop_lags=incl_subop_lags;
conf.balance_task_epochs=balance_task_epochs;
conf.lost_frac=lost_frac;
conf.do_record=do_record;
conf.video_quality=video_quality;
conf.gauss=is_gauss_kernel;
end
