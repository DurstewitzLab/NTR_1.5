
function class_multi_stats_report(order,cum_multi_stats_nonlin,multi_tests,tr_multi_tests,cum_multi_stats_lin)
%Private, auxiliar file, displaying a class-by-class statistical report for multiple
%discriminant analyses. Display will depend wether or not the "linear"
%results are provided for comparison
%
%Figure displaying % misccla. and % of divergent trajec
disp('    Displaying % missclassified, % divergent trajec per task-epoch...')
if nargin>4
    class_comp_1=[mean(cum_multi_stats_lin(:,:,1));mean(cum_multi_stats_nonlin(:,:,1))];
    err_class_comp_1=[std(cum_multi_stats_lin(:,:,1));std(cum_multi_stats_nonlin(:,:,1))]./sqrt(length(cum_multi_stats_nonlin(:,1,1)));
    %Div. trajec. (averaged over x-valid. realisations)
    class_comp_2=[mean(cum_multi_stats_lin(:,:,2));mean(cum_multi_stats_nonlin(:,:,2))];
    err_class_comp_2=[std(cum_multi_stats_lin(:,:,2));std(cum_multi_stats_nonlin(:,:,2))]./sqrt(length(cum_multi_stats_nonlin(:,1,2)));
    display_class_results(class_comp_1,err_class_comp_1,class_comp_2,err_class_comp_2,...
        order,'Multi-class discr.'),
else
    class_comp_1=mean(cum_multi_stats_nonlin(:,:,1));
    err_class_comp_1=std(cum_multi_stats_nonlin(:,:,1))./sqrt(length(cum_multi_stats_nonlin(:,1,1)));
    %Div. trajec. (averaged over x-valid. realisations)
    class_comp_2=mean(cum_multi_stats_nonlin(:,:,2));
    err_class_comp_2=std(cum_multi_stats_nonlin(:,:,2))./sqrt(length(cum_multi_stats_nonlin(:,1,2)));
    display_class_results(class_comp_1',err_class_comp_1',class_comp_2',err_class_comp_2',...
        order,'Multi-class discrim.'),
end
%
if ~all(all(isnan(multi_tests))),
    n_epochs=length(cum_multi_stats_nonlin(1,:,1));
    %Misscla.(%).
    disp('    #Iask-epochs misscla. (%). Normality tests. p(df,statistic):'),
    disp('Task-epoch           Chi2            Kolmog.-Smirnov            Lilliefords')
    for i=1:n_epochs
        disp([num2str(i),'      ',num2str(multi_tests(1,1)),'(',num2str(multi_tests(2,1)),',',num2str(multi_tests(3,1)),')',...
            '      ',num2str(multi_tests(1,2)),'(',num2str(multi_tests(2,2)),',',num2str(multi_tests(3,2)),')',...
            '      ',num2str(multi_tests(1,3)),'(',num2str(multi_tests(2,3)),',',num2str(multi_tests(3,3)),')']),
    end
    if nargin>4
        %Class-by class comparisons between O=1 and O=order
        disp('Task-epoch          T-test(O=1,O). p(df,t)          Wilcoxon rank-sum(O=1,O). p(n,ranksum)')
        for i=1:n_epochs
            disp([num2str(i),'       ',num2str(multi_tests(1,4)),'(',num2str(multi_tests(2,4)),',',num2str(multi_tests(3,4)),')',...
                '       ',num2str(multi_tests(1,5)),'(',num2str(multi_tests(2,5)),',',num2str(multi_tests(3,5)),')']),
        end
    end
    %Div.trajec. (%). Class-by class comparisons between O=1 and O=order
    %Only displayed if "full_stats" is selected.
    disp('     #Task-epochs divergent trajec. (%). Normality tests. p(df,statistic):'),
    disp('Task-epoch           Chi2            Kolmog.-Smirnov            Lilliefords')
    for i=1:n_epochs
        disp([num2str(i),'      ',num2str(tr_multi_tests(1,1)),'(',num2str(tr_multi_tests(2,1)),',',num2str(tr_multi_tests(3,1)),')',...
            '      ',num2str(tr_multi_tests(1,2)),'(',num2str(tr_multi_tests(2,2)),',',num2str(tr_multi_tests(3,2)),')',...
            '      ',num2str(tr_multi_tests(1,3)),'(',num2str(tr_multi_tests(2,3)),',',num2str(tr_multi_tests(3,3)),')']),
    end
    if nargin>4
        disp('Task-epoch          T-test(O=1,O). p(df,t)          Wilcoxon rank-sum(O=1,O). p(n,ranksum)')
        for i=1:n_epochs
            disp([num2str(i),'       ',num2str(tr_multi_tests(1,4)),'(',num2str(tr_multi_tests(2,4)),',',num2str(tr_multi_tests(3,4)),')',...
                '       ',num2str(tr_multi_tests(1,5)),'(',num2str(tr_multi_tests(2,5)),',',num2str(tr_multi_tests(3,5)),')']),
        end
    end
end


