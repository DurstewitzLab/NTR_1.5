function [r_data_ref,basic_stats,add_stats]=kfd_multiv(data_ref,labels_ref,order,reg,...
    data_val,labels_val,balance_trajec_length,full_stats_disp)
%Kernel Fisher-Discriminant (Mika et al., 00). v4.0b. 
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
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%Inputs:
%
%   data_ref (n_time x n_neurons) =  Must contain neural responses values over time.
%
%   labels_ref (n_time x 1) = Labels of task-epochs (Behavioral events) labelling the time
%                        windows in which such events occur. It is
%                        assumed that both columns contain all consecutive natural numbers
%                        until its maximum value.
%   order (1x1) = "1" for a regularized naive fisher discriminant (default)
%
%   reg  (1x1) = Regularization penalty, in units of % of the mean dot product of any two vectors.
%               (i.e. mean(mean(Kernel matrix)). If not specified, it will be automatically adjusted by leave-one-out
%               cross-validation. If order = "1", It penalizes directly the pooled covariance matrix.
%               Default value is 0.15. NOTE version 4.0b does not implement
%               the automatical regularization
%
%   data_val (n_time2 x n_dimensions (e.g. neurons)) = Same for the validation data which
%                                     will be then projected into the reference
%                                     eigenvectors (please see comments in
%                                     the code)
%
%   labels_val (n_time2 x 1) = Same for the validation data.
%
%
%   balance_trajec_length (1x1) = If >0, the length of any class-specific
%                                 trajectory is upper-bounded by this value.
%                                 In addtion, if >0,prior probabilities of each class are
%                                 set to 0.5 independently of the
%                                 number of patterns. Useful when classes have very
%                                 different sizes (not used by default).
%
%   full_stats_disp (1x1) = If >0 Computes and displays in the command window
%                      other statistics of validation set
%                      classification accuracy (besides the posterior probabilities and
%                      percentage of errors) as well as different
%                      parametric and non-parametric tests.
%                      Warning: requires Matlab "Statistics Toolbox".
%Outputs:
%    basic_stats (1x12) = Structure with the following fields
%                         IMPORTANT NOTE: All they refer to the validation
%                         set, they are out-of-sample statistics.
%
%                             basic_stats.percent_missca (1x1) = Total percent
%                               of  misclassified vectors across all classes.
%
%                             basic_stats.error (1x1)= Total outliers (%),
%                               i.e. low-dim. vectors projected in the optimal
%                               subspace whose z-score when projected exceed a
%                               certain threshold (e.g. 99.5% of the cumulative gaussian
%                               distribution). They define error bars
%                               for the % of missclasification. Please see
%                               "class_trajec.m" for more details.
%
%                             basic_stats.percent_div_trajec (1x1) = Total percent of
%                               divergent trajectories. Please see
%                               "class_trajec.m" code comments for more details.
%
%                             basic_stats.error_trajec (1x1)= Same as ".error" field but for
%                               "outlier trajectories" (in %). An "outiler trajectory" contains one or more
%                               outlier vectors. Please see below in this file comments for more details.
%
%                             basic_stats.cum_misscla_class (1 x number of classes) = Number of missclasified vectors per
%                                                                                     class.
%
%                             basic_stats.n_bins_per_class (1 x number of classes) = Number of vectors per class
%
%                             basic_stats.cum_diverg_trajec_class (1 x number of classes) = Number of divergent trajectories per class
%
%                             basic_stats.cum_trajec_class (1 x number of classes)= Number of trajectories per class
%
%                             basic_stats.predic_class (1x number of patterns (=length(labels)) = Predicted class for each pattern.
%
%                             basic_stats.shortest_mah_vec (1x number of patterns (=length(labels)) = Mahalanobis distance corresponding to the predicted class
%                                                                                                    over time
%
%                             basic_stats.posterior_probs_vectors (number of classes x number of patterns (=length(labels)) =
%                                                                                     posterior for each pattern and class (i.e
%                                                                                     "probability of classification into such
%                                                                                      class given each time pattern")
%
%                             basic_stats.lik_vectors ( max(labels) x length(labels)) = Likelihood of each pattern
%                                 given a certain class. NOTE:
%                                 not normalized across different classes for
%                                 each time bin. Please see "class_trajec.m" for more details.
%
%   add_stats (1x12) =  Warning: requires Matlab "Statistics Toolbox".
%                       IMPORTANT NOTE: All they refer to the validation
%                       set, they are out-of-sample statistics.
%                       These will be only reported and displayed in the command window if
%                       full_stats>0. Only few of these test are returned.
%                       This can be easily changed within "class_trajec.m". Please see
%                       this file and Matlab Statistical Toolbox hep
%                       for details.
%                          add_stats.Lilliefords (number of classes x 2) First column contains hypothesis verification
%                                   (h=0 means rejected) and second column
%                                   significance. Small values indicate
%                                   that distribution departure form
%                                   gaussian and thus nonparametric tests
%                                   are more advisable.
%                                   See other normality test in "class_trajec.m"
%                            Details of the following statics (plus others not returned)
%                            can be found in "class_trajec.m" code.
%
%                             add_stats.P_Manova (1 x number of classes)
%                             add_stats.Estimated_dim_Manova (1x1)
%                             add_stats.Between_class_Mahala_Manova (1 x number of classes)
%                             add_stats.Multicomp_Anova
%                             add_stats.Multicomp_Kruskal
%                             add_stats.Anova_stats
%                             add_stats.Anova_table
%                             add_stats.Kruskall_stats
%                             add_stats.Kruskall_table
%______________________________________________________________________

%I-PRE-PROCESSING. Checkings and pre-processings.
if nargin<8, full_stats_disp=0; end,
if nargin<7, balance_trajec_length=30; chance_priors=1; else chance_priors=1; end;
if nargin<5, data_val=data_ref; labels_val=labels_ref; warning('kfd_multiv:noValitation','No validation data specified! computing in-sample classification statistics'), end
if (nargin<4)||(reg<0), reg=1e-3; warning('kfd_multiv:noReg','No valid regularization specified. Set to 1e-3 of the mean kernel matrix'), end
if (nargin<3)||(order<1)||(order>10), warning('kfd_multiv:noOrder','Polynomial order not specified or invalid, set to 1 (a regularized Fisher Discriminant)'); order=1; end
if nargin<2,
    error('At least, the data matrix and task-epoch labels have to be specified'),
end
if full_stats_disp, disp(['   Discriminant analyses for order ',num2str(order),'...'] ); end
%Control over regularization
if (reg<0)
    warning('kfd_multi:regulOutOfRange',('Regularization constant not set or is out of range, set to to 1e-5 x <kernel matrix>')),
    reg=1e-3;
end
%Checking out that all labels are represented both in reference and validation sets
m_labels_ref=max(labels_ref);min_labels_ref=min(labels_ref(labels_ref>-1));l_labels_ref=length(labels_ref);
set_labs_ref=[];
for i=min_labels_ref:m_labels_ref
    if ~(all(i.*ones(l_labels_ref,1)-labels_ref))
        set_labs_ref=[set_labs_ref,i];
    end
end
n_labels_ref=length(set_labs_ref);
%
m_labels_val=max(labels_val);min_labels_val=min(labels_val(labels_val>-1));l_labels_val=length(labels_val);
set_labs_val=[];
for i=min_labels_val:m_labels_val
    if ~(all(i.*ones(l_labels_val,1)-labels_val))
        set_labs_val=[set_labs_val,i];
    end
end
n_labels_val=length(set_labs_val);
%
set_labs=set_labs_ref;
if (n_labels_val>1)&&(n_labels_ref>1)
    %Removing labels which are not in both reference and validation set
    for i=1:n_labels_ref,
        if ~any(labels_val==set_labs_ref(i)),
            warning('kfd_multi:LabelMiss',['Label ',num2str(set_labs_ref(i)),' missing in validation set']),
            %disp('Chance predictions when compared with this class will be returned');
            labels_ref=labels_ref(labels_ref~=set_labs_ref(i));
            data_ref=data_ref(labels_ref~=set_labs_ref(i),:);
            set_labs=set_labs(set_labs~=set_labs_ref(i));
        end
    end
    for i=1:n_labels_val,
        if ~any(labels_ref==set_labs_val(i)),
            warning('kfd_multi:LabelMiss',['Label ',num2str(set_labs_val(i)),' missing in reference set']),
            %disp('Chance predictions when compared with this class will be returned');
            labels_val=labels_val(labels_val~=set_labs_val(i));
            data_val=data_val(labels_val~=set_labs_val(i),:);
            set_labs=set_labs(set_labs~=set_labs_val(i));
        end
    end
    l_labels_ref=length(labels_ref);
    l_labels_val=length(labels_val);
    %m_labels=max(labels_ref);
    %min_labels=min(labels_ref(labels_ref>-1));
    n_labels=length(set_labs);
    %
    %
    %Labels and data vectors must have same number of rows.
    n_patt_ref=length(data_ref(:,1)); n_patt_val=length(data_val(:,1));
    if ~(l_labels_ref==n_patt_ref), error('Time bins of reference set labels and data differ'); end
    if ~(l_labels_val==n_patt_val), error('Time bins of validation oset labels and data differ'); end
    %
    %Now, data from the different task-epochs will be arranged on a single matrix,
    %one task-epoch block afther the other. This is necessary for the KFD algortihm
    %implementation shown in short.
    %It is however convenient to record the temporal order for showing later on
    %the flow.
    %
    n_bins_per_class_ref=zeros(1,n_labels);%Number of patterns per each class.
    %                                       That will be used later on in the
    %                                       algorithm for centering the means.
    n_bins_per_class_val=zeros(1,n_labels);
    time_order_ref=(1:n_patt_ref)';time_order_val=(1:n_patt_val)';%Contains the "real" time
    t_ref=[]; t_val=[];%Will contain time indexes, ordered per task-epochs
    cum_data_ref=[]; cum_data_val=[];%Will contain data, task-epochs ordered
    %                                (note: no need to allocate space because this loop is very short)
    class_lab_ref=[];%Auxiliary vector, containing again the task-epoch labels but where all vectors are ordered.
    class_lab_val=[];
    for i=1:n_labels
        ref_data_i=data_ref(labels_ref==set_labs(i),:);
        val_data_i=data_val(labels_val==set_labs(i),:);
        n_bins_per_class_ref(i)=length(data_ref(labels_ref==set_labs(i),1));
        n_bins_per_class_val(i)=length(data_val(labels_val==set_labs(i),1));
        cum_data_ref=[cum_data_ref;ref_data_i];%Data is acumulated in class of task-epochs
        cum_data_val=[cum_data_val;val_data_i];
        t_ref=[t_ref;time_order_ref(labels_ref==set_labs(i))];
        t_val=[t_val;time_order_val(labels_val==set_labs(i))];
        class_lab_ref=[class_lab_ref;(set_labs(i)).*ones(n_bins_per_class_ref(i),1)];
        class_lab_val=[class_lab_val;(set_labs(i)).*ones(n_bins_per_class_val(i),1)];
    end
    %Number of 'labeled' patterns will be the only ones used in the
    %following
    n_patt_ref=sum(n_bins_per_class_ref);n_patt_val=sum(n_bins_per_class_val);
    %
    %II-REFERENCE SET DISCRIMINANT. Reference set will be used to determine the best
    %discriminating directions.
    %Note: for order 1, kernels are not necessary, and a
    %simple FDA making explicit use of the covatiance matrixes can be used.
    %Kernel formulation yields to the same result.
    %
    %II.1) Compute 'gram' (kernel, K) Matrix i.e. dot product of any two high-dim
    %vectors.
    K=zeros(n_patt_ref,n_patt_ref);
    for i=1:n_patt_ref,
        for j=i:n_patt_ref,
            K(i,j)=multinom_kernel(cum_data_ref(i,:),order,cum_data_ref(j,:));
            %K(i,j)=exponential_kernel(delta,sigma,data,data);
            K(j,i) = K(i,j); %Introducing this symmetry condition is correct
            %                because kernel function is symmetric with respect
            %                to i, j (data pattern indexes).
        end
    end
    K=K./n_patt_ref;
    %
    %II.2) Compute the 'dual' to the means difference between classes ('M' matrix) and
    %of pooled covariance ('N' matrix) on a high- dimensional space.
    %
    %First of all we need a 'dual' version of the averages per class.
    %For that means, an auxiliary column of ones and zeros has to be
    %firstly constructed per class.
    hig_dim_means=cell(1,n_labels);
    cum_patterns=0;cum_means=zeros(n_patt_ref,1);
    for i=1:n_labels
        x=zeros(n_patt_ref,1);
        %Building up an auxiliary vector which contains "1" in the positions where the
        %the time bin (pattern) corresponds to the "ith-label" and zeros
        %elsewhere.
        x(cum_patterns+1:cum_patterns+n_bins_per_class_ref(i))=ones( n_bins_per_class_ref(i),1);
        %See e.g. Sch�lkopf & Smola 02 ch 14 for justification of the next step.
        hig_dim_means{i}=(1/n_bins_per_class_ref(i))*(K*x);
        cum_patterns=cum_patterns+n_bins_per_class_ref(i);
        cum_means=cum_means+hig_dim_means{i};
    end
    centroids_mean=cum_means./n_labels;%Just the mean of the class-means (or centroids)
    %
    %Again, see e.g. Sch�lkopf & Smola 02 ch 14 for justification of the next
    %step. M will represent the "dual" mean-differences (numerator in the FD ration)
    %and N the "dual" pooled covariance.
    M=zeros(n_patt_ref,n_patt_ref); N=K*K';
    for i=1:n_labels
        M=M+(hig_dim_means{i}- centroids_mean)*(hig_dim_means{i}- centroids_mean)';
        N=N-(n_bins_per_class_ref(i)*(hig_dim_means{i}*hig_dim_means{i}'));
    end
    M=M.*(n_labels-1); %across-class unbiased covariance
    %
    %Regularizing
    mK=abs(mean(mean(K)));
    if full_stats_disp>2, disp(['     Mean K is ',num2str(mK)]), end
    C=reg*mK;%The value of 'C' is such that the maximum eigenvector
    %         is proportional to alpha. This is taken as a stable solution.
    %         This value seems ok. C cannot be much smaller than
    %         0.01*mean(mean(K))in MSUA-like data. The bigger, the worse
    %         in-sample classification.
    N=N+C*K;%Emulates SVM-like penalization (i.e. complexity).
    %
    %Extracting eigenvalues and eigenvectors
    [evecs,evals]=eig(M,N);
    evals = real(diag(evals));
    %Warning: eigenvalues have to be ordered.
    [~,index]=sort(abs(evals),'descend');
    evecs=evecs(:,index);
    %Normalizing eigenvectors in the high-dim space. It is not compelling.
    %evals=evals(index(end:-1:1));
%    for i=1:n_patt_ref,
%         evecs(:,i) = evecs(:,i)/(sqrt(abs(evals(i))));
%     end
    alpha=evecs(:,1);%Is the optimum "dual" of the discriminating direction.
    %
    %II.3) Projecting data into the dual of the maximum discriminating
    %directions for plotting.
    %Gram-Schmith orthogonalization, unncoment next lines to activate.
    %This provides a more realistic plot but then axes are no longer
    %canonical variates
    v2=evecs(:,2);  v3=evecs(:,3);
    axe2=v2;%-( (alpha'*v2)/(v2'*v2) )*alpha;
    axe3=v3;%-( (alpha'*v3)/(v3'*v3) )*alpha-( (v2'*v3)/(v3'*v3) )*v2;
    r_data_ref=K*[alpha,axe2,axe3];
    %Re-constructing the real time order .
    %Note that, just for re-order, the 'unlabeled' samples must not
    %be counted, then a new 'real time' vectors of indexes has to be made
    %as follows:
    [~,used_index]=sort(t_ref);
    r_data_ref=r_data_ref(used_index,:);
    %
    %III-VALIDATION. Validation set vectors will be expanded in terms of the reference set
    %and used for classifying the validation set vectors
    %
    %III.1) Recompute now the covariane in the feature space
    K2=zeros(n_patt_val,n_patt_ref);%Note: is not a proper gram matrix anymore, is not symmetric.
    if nargin>4
        for i=1:n_patt_val,
            for j=1:n_patt_ref,%Note: no mercer kernel symmetry condition introduced here
                %              (i.e K2(i,j) <> K2(j,i), in general. Thus,
                %               "j" index starts with "1" and not with "i").
                K2(i,j)=multinom_kernel(cum_data_ref(j,:),order,cum_data_val(i,:));
                %K2(i,j)=exponential_kernel(delta,sigma,data,data);
            end
        end
        K2=K2./n_patt_val;
    else
        K2=K;
    end
    r_data_val= K2*evecs(:,1:n_labels-1);%This is the projection into the optimum discrimina-
    %                                     ting subpace as defined by the first mlabels-1
    %                                     eigenvectors. It will be used for classification
    %                                     of this validation set vectors.
    %                                     If classification is performed in
    %                                     here on a pairwise based, only the maximum
    %                                     discriminating direction 'alpha' will be used.
    %III.2) Statistical analyses of validation set
    %Recover the real temporal order
    [~,used_index]=sort(t_val);t_val=t_val(used_index);
    r_data_val=r_data_val(used_index,:);
    class_lab_val=class_lab_val(used_index);
    [basic_stats,add_stats]=class_trajec(r_data_val,class_lab_val,t_val,chance_priors,balance_trajec_length,full_stats_disp);
    %
    %If some label was not present in validation set, chance predictions
    %will be added to the real ones
    basic_stats=fill_predic(basic_stats,set_labs_ref,set_labs);
else
    warning('kfd_multiv:Single_class','Only one class is present!. No analyses returned')
    %Obviously, only perform comparisions if there are at least two labels
    basic_stats=[];
    add_stats=[];
    r_data_ref=[];
end
if full_stats_disp>1,
    disp(['   FDA over high-dimensional O=',num2str(order),' spaces finished'])
    disp('___  ___  ___'),
end

