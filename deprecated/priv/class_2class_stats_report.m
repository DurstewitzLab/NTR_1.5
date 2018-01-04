
function class_2class_stats_report(order,cum_pair_stats_nonlin,pair_tests,tr_pair_tests,combinations,cum_pair_stats_lin)
%april 011
%Private, auxiliar file, displaying all pairwise-class comparions
%statistical report for two-class discriminant analyses. Display will depend wether or not the "linear"
%results are provided for comparison
%
%Figure displaying % misccla. and % of divergent trajec.
disp('    Displaying % missclassified, % divergent trajec per task-epoch pair...')
if nargin>5
    class_comp_1=[mean(cum_pair_stats_lin(:,:,1));mean(cum_pair_stats_nonlin(:,:,1))];
    err_class_comp_1=[std(cum_pair_stats_lin(:,:,1));std(cum_pair_stats_nonlin(:,:,1))]./sqrt(length(cum_pair_stats_nonlin(:,1,1)));
    %Div. trajec. (averaged over x-valid. realisations)
    class_comp_2=[mean(cum_pair_stats_lin(:,:,2));mean(cum_pair_stats_nonlin(:,:,2))];
    err_class_comp_2=[std(cum_pair_stats_lin(:,:,2));std(cum_pair_stats_nonlin(:,:,2))]./sqrt(length(cum_pair_stats_nonlin(1,2,:)));
    display_class_results(class_comp_1,err_class_comp_1,class_comp_2,err_class_comp_2,...
        order,'pair-class discr.',combinations),
else
    class_comp_1=mean(cum_pair_stats_nonlin(:,:,1));
    err_class_comp_1=std(cum_pair_stats_nonlin(:,:,1))./sqrt(length(cum_pair_stats_nonlin(:,1,1)));
    %Div. trajec. (averaged over x-valid. realisations)
    class_comp_2=mean(cum_pair_stats_nonlin(:,:,3));
    err_class_comp_2=std(cum_pair_stats_nonlin(:,:,3))./sqrt(length(cum_pair_stats_nonlin(:,1,3)));
    display_class_results(class_comp_1',err_class_comp_1',class_comp_2',err_class_comp_2',...
        order,'Two-class discrim.',combinations),
end
if ~all(all(isnan(pair_tests))),
    %Misscla.(%)
    disp('    #Pairs of task-epochs comparisons misscla. (%). Normality tests p(df,statistic):'),
    disp('Task-epoch           Chi2            Kolmog.-Smirnov            Lilliefords')
    n_epochs=length(cum_pair_stats_nonlin(1,:,1));
    for i=1:n_epochs
        disp([num2str(i),'      ',num2str(pair_tests(1,1)),'(',num2str(pair_tests(2,1)),',',num2str(pair_tests(3,1)),')',...
            '      ',num2str(pair_tests(1,2)),'(',num2str(pair_tests(2,2)),',',num2str(pair_tests(3,2)),')',...
            '      ',num2str(pair_tests(1,3)),'(',num2str(pair_tests(2,3)),',',num2str(pair_tests(3,3)),')']),
    end
    if nargin>5
        %Class-by class comparisons between O=1 and O=order
        disp('Task-epoch          T-test(O=1,O). p(df,t)          Wilcoxon rank-sum(O=1,O). p(n,ranksum)')
        for i=1:n_epochs
            disp([num2str(i),'       ',num2str(pair_tests(1,4)),'(',num2str(pair_tests(2,4)),',',num2str(pair_tests(3,4)),')',...
                '       ',num2str(pair_tests(1,5)),'(',num2str(pair_tests(2,5)),',',num2str(pair_tests(3,5)),')']),
        end
    end
    %Div.trajec. (%)
    disp('     #Pairs of task-epochs comparisons divergent trajec. (%). Normality tests p(df,statistic):'),
    disp('Task-epoch           Chi2            Kolmog.-Smirnov            Lilliefords')
    for i=1:n_epochs
        disp([num2str(i),'      ',num2str(tr_pair_tests(1,1)),'(',num2str(tr_pair_tests(2,1)),',',num2str(tr_pair_tests(3,1)),')',...
            '      ',num2str(tr_pair_tests(1,2)),'(',num2str(tr_pair_tests(2,2)),',',num2str(tr_pair_tests(3,2)),')',...
            '      ',num2str(tr_pair_tests(1,3)),'(',num2str(tr_pair_tests(2,3)),',',num2str(tr_pair_tests(3,3)),')']),
    end
    if nargin>5,
        %Class-by class comparisons between O=1 and O=order
        %Only displayed if "full_stats" is selected.
        disp('Task-epoch          T-test(O=1,O). p(df,t)          Wilcoxon rank-sum(O=1,O). p(n,ranksum)')
        for i=1:n_epochs
            disp([num2str(i),'       ',num2str(tr_pair_tests(1,4)),'(',num2str(tr_pair_tests(2,4)),',',num2str(tr_pair_tests(3,4)),')',...
                '       ',num2str(tr_pair_tests(1,5)),'(',num2str(tr_pair_tests(2,5)),',',num2str(tr_pair_tests(3,5)),')']),
        end
    end
end

