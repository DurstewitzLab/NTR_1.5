function [basic_stats,add_stats]=class_trajec(data,labels,time,chance_priors,balance_trajec_length,full_stats_disp)
%V4.0 beta, not guaranteed.
%Warning: requires Matlab "Statistics Toolbox".
%Auxiliary file. Reports a series of statistics in a validation set:
%missclasifed vectors, divergent class-specific-trajectories etc. etc.
%
%Based on: Balaguer-Ballester E, Lapish C, Seamans JK and Durstewitz D 2011
% "Attracting Dynamics of Frontal Cortex Populations during Memory Guided Decision-Making".
% PLoS Comput. Biol., In press.
%by BCCN Heidelberg-Mannheim.
%
%Updates will be provided in this address:
%http://www.bccn-heidelberg-mannheim.de/
%Contact adresses: emili.balaguer@zi-mannheim.de
%                  daniel.durstewitz@zi-mannheim.de
%                  emili.balaguer@uv.es
%Log of changes: Please specify here your code modifications
%(#Line(s)@author: Change description). This will be useful if assistance is
%required
%#469, 508, 543@ebb:mean J-S divergence across pairs is now a new field of basic_stats
%#537-542@ebb: Catching error of Manova due to data inconsistencies 
%#478, 497@ebb: Catching test exceptions
%#571-577@ebb. Returning Wilik's lambda of the manova, as well as averaged
%                   mahalanobis distances as new fileds of "basic_stats"
%_________________________________________________________________________
%Inputs:
%
%   data (n_time x subspace dimensionality) = Neurla responses over time, spanned in an
%                                             relevant subspace (e.g. the
%                                             one from multiple
%                                             discriminant analysis).
%
%   labels (n_time x 1) = Labels of Task-epochs (Behavioral events) labelling the time
%                        windows in which such events occur. It is
%                        assumed that both columns contain all number from 1 to
%                        its maximum value.
%
%   chance_priors (1x1) = If >0, prior probabilities of each class are 0.5 independently of the
%                         number of patterns. Useful when classes have very
%                         different sizes.
%
%   balance_trajec_length (1x1) = If >0, the length of any class-specific
%                                 trajectory is upper-bounded by this value.
%                                 Useful when classes have very
%                                 different sizes (not used by default).
%
%   full_stats_disp (1x1) = If >0 Computes and displays in the command window
%                      other statistics of validation set
%                      classification accuracy (besides the posterior probabilites and
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
%                           See code for more fields returned on
%                           basic_stats (Mean Jensen Shannon, summed WIlk's
%                           lambda, mean Mahalanobis between means in a
%                           MANOVA analysis)
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
%__________________________________________________________________________

%I. PRE-PROCESSING. Checkings and data arrangement
%
%I.1 Checkings
if nargin<4, balance_trajec_length=0; end %For much longer trajectories of one epoch with respect to the others, one migth want to cut them in segments of equal duration
if nargin<3, error('Validation dataset, the corresponding class-labels and the real time column have to be specified'), end
[n_patt,subspace_dim]=size(data);
m_labels=max(labels);min_labels=min(labels(labels>-1));l_labels=length(labels);
if ~(l_labels==n_patt), error('Class labels and dataset times must have the same length'), end
set_labs=[];
for i=min_labels:m_labels
    if ~(all(i.*ones(l_labels,1)-labels))
        set_labs=[set_labs,i];
    end
end
%Current warning state
warn_state=warning('query','all');
n_labels=length(set_labs);
%
%I.2 Hard-codded settings
no_null_resp=6;%Minimum number of no-null time bins per each class.
confidence_class=99.8;%In percent: threshold for detecting if a concrete pattern is an outliner,
%i.e. over aprox 3.09 z-scores has an ideal probability to occur of 0.001 (99.8% conf); or 2.576 for a
%prob 0.005 (99.5%) assuming gaussianity of the projected data (see those numbers in any statistical
%book e.g. Makridakis et al., 1999). This will serve to provide error bars different
%to the typical SEM i.e. without making any average
z_lim=norminv([0 0.01*confidence_class],0,1);z_lim=z_lim(2);%positive z-score for the specified confidence
%
min_trajec_length=2;%Minimum length of a trajectory (in time bins) in order to be considered
%                   as either convergent or divergent when not all vectors are classified in
%                   whithin the same (either correctly or incorrectly)
%                   class. Please see next section for more details.
%
gauss_factor=(2*pi)^(-subspace_dim/2);%The probabilites of classification will be
%                                       computed in an optimum
%                                       discriminating subpace, which
%                                       can be of any dimension. For
%                                       instance, if it was obtained
%                                       via Fisher discrimnant then the
%                                       dimensionality of such subspace
%                                       will be number of classes minus
%                                       one (i.e. "m_labels-1").
%                                       This is the constant factor
%                                       of any multivatriate normal
%                                       distribution.
%
%I.3 Data arragement and basic statistics
%Split data and compute basic statistics separately for each class
%Note: "cell" variables used due to the multi-dimensionality of data.
n_bins_per_class=zeros(1,n_labels);
mean_l_trajec=zeros(1,n_labels);
std_l_trajec=zeros(1,n_labels);
class_data=cell(1,n_labels);
class_times=cell(1,n_labels);
priors=0.5.*ones(1,n_labels);means=cell(1,n_labels);covariances=cell(1,n_labels);
%Note: next loop is very small thus does not slows down significantly the
%code
for i=1:n_labels,
    %Split data into classes
    class=set_labs(i);
    class_data{i}=data(labels==class,:);
    class_times{i}=time(labels==class);
    n_bins_per_class(i)=length(class_times{i});
    %Now control the fact that we have any non-almost-null vector
    if length(  find(    (abs(class_data{i}))>0  )    )<no_null_resp,
        warning('class_trajec:nullResp',['Unreliable analyses: Less than ',num2str(no_null_resp),' no-null responses'])
    end
    %Set priors, means and covariance matrix per class (note the
    %multi-dimensionality of such quantities). Equal priors are convenient when
    %there are many more samples in one class whith respec to the others.
    if chance_priors>0,
        priors(i)=0.5;
    else
        priors(i)=n_bins_per_class(i)/n_patt;
        warning('class_trajec:chancePriors','Chance priors not set'),
    end
    means{i}=mean(class_data{i});
    covariances{i}=cov(class_data{i});
    %Compute the mean duration of trajectories per class. That will be useful
    %in the next section.
    trajec_pointer=0;trajec_lengths=[];%Upper bound in the number of trajectories
    n_trajects=0;                                                %within the ith class. For saving
    %                                                            memory (next loop can be computationally
    %                                                            demanding).
    current_class_times=class_times{i};
    for j=1:n_bins_per_class(i),
        %Are we within the same class-specific-trajectory?
        
        if (j==1)||(current_class_times(j)==1+current_class_times(j-1))
            trajec_pointer=trajec_pointer+1;
        else
            n_trajects=n_trajects+1;
            trajec_lengths=[trajec_lengths,trajec_pointer];
            trajec_pointer=1;
        end
    end
    n_trajects=n_trajects+1;
    trajec_lengths=[trajec_lengths,trajec_pointer];
    trajec_pointer=1;
    mean_l_trajec(i)=mean(trajec_lengths(trajec_lengths>0));%Note that zero vectors are not considered
    %sem_l_trajec(class)=std(trajec_lengths(trajec_lengths>0))/sqrt(n_trajects);
    std_l_trajec(i)=std(trajec_lengths(trajec_lengths>0));
    if full_stats_disp>1
        disp(['     Mean trajectory length for class ',num2str(class),...
            '=',num2str(mean_l_trajec(i)),', SEM=',...
            num2str(std_l_trajec(i)./sqrt(n_trajects)),' time bins']),
    end
end
%
%It could be one trajectory is much longer than the others. In this
%situation, balance_trajec_length>0 causes that trajects. much longer than
%the others are cut-off.
max_trajec_length=inf;
if balance_trajec_length>0,
    max_trajec_length=max(mean_l_trajec);
    %max_trajec_length=m_mean+(balance_trajec_length*mean(std_l_trajec));
    if max_trajec_length>balance_trajec_length, max_trajec_length=balance_trajec_length; end
    if full_stats_disp>1, disp(['     Maximum trajectory length restricted to ',num2str(max_trajec_length),' time bins']), end
end
%
%Finally, this auxiliar operation (filling 'class data' with NaN so that all classes have the same length)
%will be beneficial for further statistical testing (secion III)
max_patterns=max(n_bins_per_class);
for i=1:n_labels,
    if ~(max_patterns==n_bins_per_class(i)),
        class_data{i}=[class_data{i};NaN.*ones(max_patterns-n_bins_per_class(i),subspace_dim)];%Note, that is the dimensionality of class_data
    end
end
%
%II-BAYES CLASSIFICATION. Naive bayesian classification and convergency of this clasification
%for each trajectory
%
%An incorrectly classified point could be due to the fact that bayesian criterion of the
%discriminat is nor realistic enough for segregating several attractors (in an out-of-sample validation sense of course)by a linear decision rule.
%However, they migth be still separable by an unknown boundary.
%Convergency test is as follows: a wrongly classifed point does belong to a convergent
%trajectory only if all of the future points for now to the end of the
%trajectory will become well-classified.
%Otherwise, the entire trajectory will be considered as divergent.
%
misscla=0;%Total number of missclasified vectors
misscla_consec=0;%This variable measures the total number of incorrect points which has occurred in sequence within the
%                  same trajectory
outliers=0;%Total number of outliers (see section I)
outlier_trajecs=0;
diverg_trajec=0;%Total number of divergent trajectories (see upper pargraph)
already_outlier=0; %flag, auxiliar ariable to label a trajectory as outlier
trajec_pointer=1;%Points to the position in the trajectory (current length)
cum_lik_mah=zeros(n_labels,n_patt);%Columns contain mahalanobis distances of each pattern from
%                                   the centroid of the "true" class
%                                   where it belongs.
cum_lik_probs=zeros(n_labels,n_patt);%Columns contain the likelihood of
%                                     each ith pattern given a certain
%                                     class (i.e. patterns in columns).
%                                     Normalization is such that sum across
%                                     all infinite ith-patterns of a
%                                     class "c"=1. See below
cum_posterior_probs=zeros(n_labels,n_patt);%Columns contain the posterior probability of
%                                           a certain class to occur given each ith
%                                           pattern (i.e. patterns in columns).
%                                           Warning: each column is
%                                           not normalized such that
%                                           this pattern has a
%                                           probability "1" of been
%                                           classified into one class.
%                                           See below
cum_predic_class=zeros(1,n_patt);%More likely class for each time pattern
cum_max_logs=zeros(1,n_patt);%log-likelihood for each time pattern
cum_misscla_class=zeros(1,n_labels);
cum_diverg_trajec_class=zeros(1,n_labels);
cum_trajec_class=zeros(1,n_labels);
n_trajects=0;
%
for i=1:n_patt,
    c=labels(i);%This is the experimentally correct class for ith-pattern
    %Are we in the same class-specific-trajectory?
    if (i==1)||((c==labels(i-1))&&(trajec_pointer<=max_trajec_length))
        trajec_pointer=trajec_pointer+1;
    else
        %Convergency analysis of the just-finished trajectory:
        %If we have "jumped" to the next trajetory then check out
        %if the last (or lasts) vectors of the just-ended previous
        %trajectory where incorrectly classfied.
        %If that is the case, then the trajectory can be potentially
        %"diverging" form the current "class" as will be label as such.
        %
        %Refinement: Class-specific trajectories made of less than
        %"min_trajec_length" vectors (e.g. those made of just two time bins)
        %cannot really be considered as divergent unless all those vectors where
        %incorrectly assigned
        if (i>1)&&...
                (( misscla_consec>0)&&(trajec_pointer>min_trajec_length))%||( misscla_consec==trajec_pointer)
            diverg_trajec=diverg_trajec+1;
            cum_diverg_trajec_class(set_labs==c)=cum_diverg_trajec_class(set_labs==c)+1;
        end
        %Too short class-specific trajectories are not convergent either;
        %then, they should not be reported anymore (unless all of their points become
        %correctly classified, as in with divergent trajectories).
        if (trajec_pointer>min_trajec_length)%||(misscla_consec==0)
            n_trajects=n_trajects+1;%Re-counting the number of trajectories per class
            cum_trajec_class(set_labs==labels(i-1))=cum_trajec_class(set_labs==labels(i-1))+1;
        end
        trajec_pointer=1;%Trajectory length of the newly started trajectory
        misscla_consec=0;% In this brand new trajectory there are neither incorrectly
        %                 or correctly classified points yet
        already_outlier=0;%This will label the "outliner" trajectories
    end
    %
    %Compute log of the P(class "j"|for each firing pattern)
    %for each "jth-class" (assuming gaussianity; therefore Bayes criterion becomes optimal).
    % This is probably computationally unneficient (algorihm has now
    %(n_patterns*(m_labels^2)) iterations) but the implementation is
    %clearer (please contact emili.balaguer@zi-mannheim.de for a faster implementation).
    joint_prob=zeros(1,n_labels);lik_prob=zeros(1,n_labels);z_square=zeros(1,n_labels);
    log_score=zeros(1,n_labels);
    for j=1:n_labels
        %Recovering second-order statistics for each class
        %c_n=n_bins_per_class(class);
        j_prior=priors(j);
        %j_mean_l_trajec=mean_l_trajec(j);
        j_mean=means{j};
        j_cov=covariances{j};
        %Argment of the exponential function for the ith-pattern
        %belonging to class "c" (which is c_data(i,:))
        %"as if" he would belong to the undelying gaussian
        %distribution of class "j". This is just a squared z-score.
        %
        %jcov matrix is often close to singular. Next code alleviates the
        %problem arising when inverted
        %
        regul_attempts=10; lastwarn('No warning');%Contains the warning message;
        for r=1:regul_attempts
            tol=(r-1)*exp(r-1).*eps;
            j_cov=j_cov+tol.*mean(mean(j_cov)).*eye(subspace_dim);%Avoid numerical problems. up to eps*1e-4 in Demo data.
            z_square(j)=((data(i,:)- j_mean)/j_cov)*(data(i,:)- j_mean)';
            [is_singular,id]=lastwarn;
            if (r>2)&&(full_stats_disp>1), warning('class_trajec:inverted_cov',['regularizing covariance, tol.=',num2str(tol)]), end
            %If the last warning start with "Matrix", do reffer to the
            %singularity. Then, reset warning message
            if (length(is_singular)>14)&&(strcmp(is_singular(1:14),'Matrix is clos')),
                previous=id;
                lastwarn('No warning');
            elseif(length(is_singular)>14)&&(strcmp(is_singular(1:14),'Matrix is sing')),
                warning('off',id);
            else
                if (r>1),
                    if (full_stats_disp>1), warning('class_trajec:inverted_cov',['Succesful inversion of covariance matrix in the discrim. subspace (tol=',num2str(tol),'<COV>)']); end
                    warning('off',previous), warning('off','class_trajec:inverted_cov'),
                end
                break
            end
        end
        %Now obtain the probability P(firign pattern "i" of class "c"|class "j")
        %Note that likelihoods where not normalized across groups i.e.
        %it is is not assumed that the pattern will for certain be
        %classifed.
        lik_prob(j)=gauss_factor*((det(j_cov))^-0.5)*exp(-0.5*z_square(j));
        joint_prob(j)=lik_prob(j)*j_prior;
        log_score(j)=log(joint_prob(j));
    end
    %Predicted class
    [max_log,predic_j]=max(log_score);
    predicted_class=set_labs(predic_j);
    if predicted_class==c,
        misscla_consec=0;
    else
        misscla=misscla+1;
        misscla_consec=misscla_consec+1;
        cum_misscla_class(set_labs==c)=cum_misscla_class(set_labs==c)+1;
    end
    if abs(sqrt(z_square(predic_j)))>z_lim,
        outliers=outliers+1;
        %Current trajectory will be labeled as 'outlier_trajec'
        if already_outlier==0;
            outlier_trajecs=outlier_trajecs+1;
            already_outlier=1;
        end
    end;
    %Those will serve to create error bars without any average.
    posterior_prob=joint_prob./sum(joint_prob);%True posterior
    cum_posterior_probs(:,i)=posterior_prob;
    cum_lik_probs(:,i)=lik_prob;
    cum_lik_mah(:,i)=z_square;
    cum_max_logs(i)=max_log;
    cum_predic_class(i)=predicted_class;
end
%Restoring warning state
warning(warn_state);
%
%III-STATISTICAL ANALYSES. Tests (param & nonparam.) and comparisons between pdf's of each group
%
%III.1-Basic statistics
percent_misscla=100*(misscla/n_patt);conf_int= 100*(outliers)/n_patt/2;
percent_div_trajec=100*(diverg_trajec/n_trajects);conf_int_trajec=100*(outlier_trajecs)/n_trajects/2;
if percent_misscla>(((n_labels-1)/n_labels)*100),
    %   percent_misscla=((n_labels-1)/n_labels)*100; percent_div_trajec=((n_labels-1)/n_labels)*100;
    percent_misscla=50; percent_div_trajec=50;
    warning('class_trajec:wrongClassific','Unsuccesful classification. Returning chance values'),
end
%
basic_stats.percent_misscla=percent_misscla;%missclasified (%)
basic_stats.error=conf_int;%outliers (%)
basic_stats.percent_div_trajec=percent_div_trajec;%divergent (%)
basic_stats.error_trajec=conf_int_trajec;%outlier trajectories (%)
basic_stats.cum_misscla_class=cum_misscla_class;%missclasified vectors per class (number)
basic_stats.n_bins_per_class=n_bins_per_class;%number of vectors per class
basic_stats.cum_diverg_trajec_class=cum_diverg_trajec_class;%divergent trajectories per class (number)
basic_stats.cum_trajec_class=cum_trajec_class;%number of trajectories per class
%
%Next, the more detailed vector-by-vector statistics (columns) for each class (rows),
%used for obtaining upper percentages
%
basic_stats.predic_class=cum_predic_class;
basic_stats.posterior_probs_vectors=cum_posterior_probs;
basic_stats.lik_vectors=cum_lik_probs;
basic_stats.shortest_mah_vec=cum_lik_mah;
%
%Command window display of basic statistics
if full_stats_disp>0
    disp(' ')
    disp(['Classification accuracy (',num2str(subspace_dim),' discriminating axes used)[',num2str(confidence_class),'% confidence intervals]']),
    disp(['     ',num2str(percent_misscla),'[',num2str(conf_int),'] % missclasified vectors']),
    disp(['     ',num2str(percent_div_trajec),'[',num2str(conf_int_trajec),'] % divergent trajectories']),
end
%III.2 add statistics
%
%Inisializations:
%Error bar for KL-Divergences:
%p=0.001; confidence99.8% "gaussian" error bars
p=0.005; confidence_kl=99;
%p=0.01; confidence=98%
bound=p.*ones(n_patt,1);
k_l=zeros(n_labels,n_labels);
j_s=zeros(n_labels,n_labels);j_s_error=zeros(n_labels,n_labels);
%Gaussianity tests (univariate, only defined for the first data axis)
chi2_class=zeros(n_labels,2);ks_class=zeros(n_labels,2);lill_class=zeros(n_labels,2);
%Pairwise comparison tests (univariate, only defined for the first data axis)
ttest_class=zeros(n_labels,n_labels,2);ranksum_class=zeros(n_labels,n_labels,2);
ks2_class=zeros(n_labels,n_labels,2);
%Auxiliary variable that will be used later on for global testing
data_global_testing=NaN.*ones(max_patterns*n_labels,subspace_dim);
groups=NaN.*ones(max_patterns*n_labels,1);
JSs=zeros(n_labels*(n_labels-1),1); k=1; %Number of different J-S divergences (see below)
for i=1:n_labels
    %III.2.1 Univariate tests of gaussianity (only defined for the
    %most responsive axe of the optimum suspace)
    current_data=class_data{i};
    [~,selected_dim]=max(mean(current_data(~isnan(current_data(:,1)),:)));
    %Less than two pattern per class do not allow any statistic
    if length(current_data(~isnan(current_data(:,selected_dim)),selected_dim))>4,
        %Now look at each one of subspace dimensions and test for
        %gaussianity
        try
            %Chi2 test of gaussianity
            [chi2_class(i,1),chi2_class(i,2)]=chi2gof(current_data(:,selected_dim));
            %Kolmogorov-Smirnov. Note, distributions has to be standarized
            %first.
            [ks_class(i,1),ks_class(i,2)]=kstest(current_data(:,selected_dim));
            %Lilliefors
            [lill_class(i,1),lill_class(i,2)]=lillietest(current_data(:,selected_dim));
            %
        catch
            warning('class_trajec:GaussTests_error','Data inconsistency: Tests of gaussianity not computed')
            chi2_class(i,1)=0; chi2_class(i,2)=0;
            ks_class(i,1)=0; ks_class(i,2)=0;
            lill_class(i,1)=0; lill_class(i,2)=0;
        end
    end
    for j=1:n_labels
            current_data_bis=class_data{j};
            %II.2.2 Univariate class-paiwise comparison tests
            try
                [ttest_class(i,j,1),ttest_class(i,j,2)]=ttest2(current_data(:,selected_dim),current_data_bis(:,selected_dim));
                [ranksum_class(i,j,1),ranksum_class(i,j,2)]=ranksum(current_data(:,selected_dim),current_data_bis(:,selected_dim));
                [ks2_class(i,j,1),ks2_class(i,j,2)]=kstest2(current_data(:,selected_dim),current_data_bis(:,selected_dim));
            catch
                ttest_class(i,j,1)=0; ttest_class(i,j,2)=0;
                ranksum_class(i,j,1)=0; ranksum_class(i,j,2)=0;
                ks2_class(i,j,1)=0; ks2_class(i,j,2)=0;
                warning('class_trajec:PairTests_error','Data inconsistency: Pairwise comparisons tests not computed')
            end
            pdf_i=cum_posterior_probs(i,:)'; pdf_j=cum_posterior_probs(j,:)';
            %
            %II.2.3 Kullbach-Leibler(k_l)and Jensen-Shannon (j_s) distances between pairs of groups
            %i.e. between the probabilities of belonging to such group
            %(univariate)
            M=0.5*(pdf_i+pdf_j);
            KL_i=kl_diverg(pdf_i,M);
            KL_j=kl_diverg(pdf_j,M);
            KL_up_i=kl_diverg(pdf_i+bound,abs(M-bound));
            KL_up_j=kl_diverg(pdf_j+bound,abs(M-bound));
            KL_down_i=kl_diverg((pdf_i-bound),(M+bound));
            KL_down_j=kl_diverg((pdf_j-bound),(M+bound));
            k_l(i,j)=kl_diverg(pdf_i,pdf_j);
            j_s(i,j)=0.5*(KL_i+KL_j);
            if i<j, 
                k=k+1; 
                if isnan(JSs(k))|| (abs(JSs(k))>10)
                    warning('class_trajec:NaN_divergence_error', [' Jensen-Shanon is not number please trace the comparison (i,j)=',num2str(i),',',num2str(j)]);
                    error([' Jensen-Shanon is not number please trace the comparison (i,j)=',num2str(i),',',num2str(j)]);
                else
                    JSs(k)=j_s(i,j);
                end
            end %For a posterior average of different JS's
            j_s_up=0.5*(KL_up_i+KL_up_j);j_s_down=0.5*(KL_down_i+KL_down_j);
            j_s_error(i,j)=max(abs(j_s(i,j)-j_s_up),abs(j_s(i,j)-j_s_down));
    end
    %Auxiliary step: data arrangement for global testing
    data_global_testing( (i-1)*max_patterns+1:i*max_patterns,:)=current_data;
    groups((i-1)*max_patterns+1:i*max_patterns)=i.*ones(max_patterns,1);
end
%
%II.2.4 Global tests, both uni- and multivariate
%JS storing
JSs=JSs(JSs>0);
basic_stats.js=mean(JSs);
%Univariate global statistics
try
    [anova_p,anova_table,anova_stats]=anova1(data_global_testing(:,selected_dim),groups,'off');
    [kruskal_p,kruskal_table,kruskal_stats]=kruskalwallis(data_global_testing(:,selected_dim),groups,'off');
    %Multiple comparisons. Of course, there are many different statistics to
    %be used, here we just report the default 'tukey-kramer' statistics implemented in
    %matlab stats. toolbox.
    %Please type ">help multcompare" for other options.
    Multicomp_anova = multcompare(anova_stats,'display','off');
    Multicomp_kruskal = multcompare(kruskal_stats,'display','off');
catch 
    warning('class_trajec:multiv_error','Data inconsistency: Multiple comparisons tests not computed. Cascade of errors is likely when displayed full statistics')
    anova_stats=-1;
    anova_table=-1;
    kruskal_stats=-1;
    kruskal_table=-1;
    Multicomp_anova=-1;
    Multicomp_kruskal=-1;
end
%
%Multivariate comparisons: MANOVA. Please type ">help manova1" for
%other statistics
try 
[D,P,manova_stats]=manova1(data_global_testing,groups);
catch 
    D=1; P=1; manova_stats.lambda=1; manova_stats.gmdist=inf;
    warning('class_trajec:Manova_error','Data inconsistency: Manova not computed')
end
if (D<1),  
    D=1; P=1; manova_stats.lambda=1; manova_stats.gmdist=inf; 
    warning('class_trajec:Manova_error','Data inconsistency: Manova not computed')
end
%Wilk's lambda, summed across the D-significant dimensions
wl=manova_stats.lambda; 
wl=wl(1:D);
basic_stats.lambda=sum(wl);
basic_stats.lambda2=wl(1)/D;
%Mahalanobis distances between group means, averaged
mhd=manova_stats.gmdist; n=length(mhd(:,1));
mhd=sum(sum(triu(mhd)))./(   ((n^2)+n)/2 );
basic_stats.mhd=mhd;
%
add_stats.Between_class_Mahala_Manova=manova_stats.gmdist;
add_stats.P_Manova=P;
add_stats.Estimated_dim_Manova=D;
add_stats.lambda=manova_stats.lambda;
add_stats.Multicomp_Anova=Multicomp_anova;
add_stats.Multicomp_Kruskal=Multicomp_kruskal;
add_stats.Anova_stats=anova_stats;
add_stats.Anova_table=anova_table;
add_stats.Kruskall_stats=kruskal_stats;
add_stats.Kruskall_table=kruskal_table;
add_stats.Lilliefords=lill_class;
%
if full_stats_disp>1%Bug in previous version
    %Command window display of add statistics and fill the otput argument
    disp(' ')
    disp('   Fulls statistics report:'),
    disp('   Note: Please see Matlab Stats. Toolbox help for further interpetation')
    disp('         Please see code for accesing to test statisitcs and degrees of freedom'),
    disp(' '),
    disp('    #Kullbach_Leibler divergence (class i, class j)'),disp(num2str(k_l));
    disp('    #Jensen-Shannon (J-S) diverg. (class i, class j)'),disp(num2str(j_s)),
    disp(['    #J-S errors at ',num2str(confidence_kl),'%=']), disp(num2str( j_s_error)),
    disp(' '),
    disp('    #Normality tests per class (h=0 accepts gaussianity., p=alpha sign.):'),
    disp('        Chi2(h,p)      Kolmog.-Smirnov(h,p)      Lilliefords(h,p)')
    for i=1:n_labels
        disp(['    ',num2str(set_labs(i)),'  (',num2str(chi2_class(i,1)),',',num2str(chi2_class(i,2)),')',...
            '       (',num2str(ks_class(i,1)),',',num2str(ks_class(i,2)),')',...
            '              (',num2str(lill_class(i,1)),',',num2str(lill_class(i,2)),')']),
    end
    disp(' '),
    disp('    #Univariate pairwise comparisons(p=alpha sign., if p<0.05 means are not equal):'),
    disp('     Note: Data projected into the maximum discriminating direction'),
    disp('     Unpaired t-test h='),disp(num2str(ttest_class(:,:,1))),disp('  p='),disp(num2str(ttest_class(:,:,2))),
    disp('     Wilcox-ranksum test h='),disp(num2str(ranksum_class(:,:,1))),disp('  p='),disp(num2str(ranksum_class(:,:,2))),
    disp('     Two-sample-Kolmog.-Smir. test h='),disp(num2str(ks2_class(:,:,1))),disp('  p='),disp(num2str(ks2_class(:,:,2))),
    disp(' ')
    disp('    #Univariate multiple comparisons(p=alpha sign.):'),
    disp(['    One-way Anova (p=',num2str(anova_p),') table:']),disp(anova_table(:,[3,5,6])),
    disp(['    Kruskall-Wallis (p=',num2str(kruskal_p),') table:']),disp(kruskal_table(:,[3,5,6])),
    disp(' ')
    disp('    #Post-hoc tests:')
    disp('    Multiple comparisons Anova'),disp(num2str(Multicomp_anova)),
    disp('    Multiple comparisons Kruskal-Wallis'),disp(num2str(Multicomp_kruskal)),
    disp(' ')
    disp('    #Multivariate multiple comparisons (Manova):'),
    disp(['     Estimated dimension=',num2str(D)']),
    disp('     Alpha sign. per dimension'),disp(num2str(P)),
    disp('     Wilks Lambda'),disp(num2str(manova_stats.lambda)),
    disp('     Between-classes means Mahalanob. distances'),disp(num2str(manova_stats.gmdist)),
    disp(' ')
    disp('   End of full statistics report'),
    disp(' ');
end
end

