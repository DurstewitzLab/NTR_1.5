function results=ntr(filename,conf)
%Starting function for the utility "neural activity trajectories reconstruction"
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
%NOTE: For a first approach to this toolbox please read 'readMeFirst.pdf'
%      For running a short demonstration of this toolbox please type ">ntr;"
%
% Requirements: Matlab v7.11 (2010b), Statistics Toolbox v7.4 (2010b).
% It can be easily adapted to work on earlier matlab versions (please contact authors if any difficulty).  
% The code runs either in Windows 7 or in Linux, but please remembers to change the path in the first code lines 
% of this file. 
%
%Performs time-series-based statistical analyses on multivariate
%activity recordings during a cognitive task, which contains separate
%perceptual or behavioral "epochs". 
%Therefore it is assumed that the experimentalist has labelled the task according 
%to separate cognitive epochs (e.g. an stimulus
%presentation, a successful choice, reward acquisition, movement to a
%definite place like an arm). 
%It is also assumed that all multivariate
%neural responses are time-ordered and simultaneously recorded.
%Trajectory convergence or divergence will be evaluated in those states defined by the 
%separate time-windows corresponding to different epochs of the task.
%Analyses will indicate whether those states behave like attracting regions of neural 
%activity or not. For statistical evaluation of such temporal properties 
%in extremely sparse, high-dim. neural interactions spaces kernel
% methods were needed. 
%__________________________________________________________________________
%
%I-HOW TO
%
%I.1) Drop a ".mat" or file (e.g. named "file_name.mat") on "./data" folder, 
%where "./" indicates the directory containing this entire toolbox or 
%"deployment directory" (Note: use ".\" for Linux). This "file_name.mat" 
%should contain a matlab matrix named "Data", of dimension (number
%of time bins x (dimensionality of neural responses+3) ), where "dimensionality 
%of neural responses" refers e.g. to the number of simultaneously recorded neurons,
%voxels, electrodes etc. and/or delayed version of them. 
%The dataset structure is the following:
%
%	Columns "1:end-3" of "Data" matrix:  Must contain neural responses over time.
%	Column "end-2" of "Data" matrix:  Must contain natural numbers>0, labelling 
%                                     the different stimulus or behavioural "epochs"
%                                     in which the experimentalist segments the task. 
%                                      "-1" encodes "no-labelled" time-bins.
%	Column "end-1" of "Data" matrix:  Must contain natural numbers>0, they are 
%                                     alternative labelling used only for trial 
%                                     trajectory display (number of   labels has to
%                                     be smaller than 8. See section code for more info).
%                                     Those labels typically represent "phases" of the
%                                     task, containing different "epochs". 
%                                     If all "epochs" are to be displayed, 
%                                     this column has to be a copy of "end-2" one.
%	Last column of "Data" matrix: Must contain natural numbers>0, labelling the 
%                                 different trials of the task. 
%                                 Trials typically represents repetitiesta por una combinación de algunos ions of
%                                 the experiment, containing each the same "phases" 
%                                  of the task.
%   
%
%I.2) Open "kspaces_configuration.m" and setup the configuration parameters.
%     Please find a detailed description of the parameters by typing 
%     ">help kspaces_config".
%I.3) Once in the deployed directory, type ">ntr("file_name");".
%     Alternatively, one can load in the workspace the data matrix variable
%     formatted as indicated in I.1 (e.g. "Data_matrix_name") and 
%     type ">ntr(Data_matrix_name);".
%
%
%As indicated earlier, for an introductory demonstration of this toolbox please type ">ntr;".
%
%II-OUTPUT ANALYSES.
%
%Plots and command window displays are controlled by parameters specified in "kspaces_config.m" file. 
%They are subdivided (for clarity purposes) into "Basic" and "Advanced" parameters.
%In this brief note, only few parameters are described.
%Please type ">help kspaces_config" and ">help setup_config" 
%for more information. Please also find detailed comments in each file of the toolbox.
%The output of this function is the next variable:
%
%results:(1x10) structure with the following fields:
%               multiClass_err = (1x1) Percentage of missclassified vectors using a
%                               multiple discriminant analysis, averaged
%                               across cross-validation blocks.
%               multiClass_diver =(1x1) Percentage of divergent trajectories using a
%                                 multiple discriminant analysis, averaged
%                                across cross-validation blocks.
%               twoClass_2mean = (1x1) Percentage of missclassified vectors using a
%                                two-class discriminant analysis, averaged
%                                across cross-validation blocks and across
%                                 pair-of-classes comparisons.
%               twoClass_2std_err=(1x1) Standard mean error of upper variable.
%               twoClass_2mean_div=(1x1) Percentage of divergent trajectories using a
%                                two-class discriminant analysis, averaged
%                                across cross-validation blocks and across
%                                pair-of-classes comparisons.
%                twoClass_2std_div=(1x1)Standard mean error of upper variable.
%                twoClass_mean_err=(1 x number of task-epochs) Percentage of missclassified
%                                vectors using a two-class discriminant analysis,
%                                averaged across cross-validation blocks and for
%                                each task-epoch pair comparison
%               twoClass_std_err= (1 x number of task-epochs) Standart mean error 
%                                 of upper variable.
%               twoClass_mean_div= (1 x number of task-epochs) Similar to
%                                   "twoClass_mean_err" but for % of divergent trajectories.
%               twoClass_std_div=(1 x number of task-epochs) Standard mean error 
%                                 of upper variable.
%__________________________________________________________________________

%Preliminaries
close all;if nargin<2, ts=tic; end
warning('on','all');
%Adding current directory to the path
addpath(cd,'-begin');%Contains the major functions. They can be modified by the user
addpath([cd,'/priv'],'-begin');%Contains auxilar utilities under development, not recommended to be modified
addpath([cd,'/data'],'-begin');%Include in here the datasets
%
disp('Neural trajectory reconstruction Toolbox 1_5')
disp(['By Emili Balaguer-Ballester. School of Engineering,'...
    ,' Bournemouth University']),
disp('& BCCN Heidelberg-Mannheim. Please contact eb-ballester@bournemouth.ac.uk for questions')
disp(' '),
%If no input, demo file is loaded
if nargin<1, filename='demo'; end
%Detecting if input is a file or a variable. If a variable, will be "data"
if any(isletter(filename)), load(filename); else data=filename; filename='workspace data matrix'; end
%See" kspaces_configuration.m" for information 
%about parameters. 
%Setting default configuration for demo file.
%Else, loading the user-specified configuration from
%"kspaces_configuration.m"
if strcmp(filename,'demo')>0,
    warning('ntr:loadingFile','Demo file used, default configuration loaded...'), 
    conf=setup_config(5,0.25,0,0,5,5,[4],1,0,-1,1,30,1,30,1,30,.1,0,100,0);               
    %Please see "setup_config.m" for more info.
    %
else
    disp(['Loading ',num2str(filename),'...'])
    if nargin<2,
        disp('Using configuration file "kspaces_configuration.m"....')
        conf=kspaces_config();
    end
end
%Data shuffling (producing a single bootstrap sample)
if (conf.shuff_epoch) || (conf.shuff_within_epoch)
    data=shuff_data(data,conf.shuff_epoch,conf.shuff_within_epoch);
end
%Figures and movie
if (conf.trial_disp>0)||(~isempty(conf.trials_flow_disp))
    visualization(data,conf.order,conf.trial_disp,conf.fixed_lag,conf.do_dcm,conf.make_video,...
        conf.trials_flow_disp,conf.balance_task_epochs,conf.regul,...
        conf.incl_subop_lags,conf.lost_frac,conf.do_record,conf.video_quality);
end
%Cross-validation
if (conf.multi>0)||(conf.pairwise>0),
results=kfd_cross_val(data,conf.order,conf.regul,conf.xval_type,conf.multi,conf.pairwise,...
    conf.full_stats_disp,conf.balance_trajec_length,conf.do_dcm,conf.incl_subop_lags,...
    conf.balance_task_epochs);
else
    results=[];
    disp('No analyses performed'),
end
%path(current_path);%Unncoment this and first lines of this function to
%                    recover previous path
%Display elapsed time only when external config is not provided 
if nargin<2
    minut=floor((toc(ts)/60));secon=floor(toc(ts)-minut*60);decim=round(toc(ts)-minut*60-secon);
    disp(' '),
    disp(['Finished in ',num2str(minut),' minutes ',num2str(secon),' seconds ',num2str(decim),' decims.']),
end




