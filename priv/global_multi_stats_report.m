
function global_multi_stats_report(order,global_multi_stats_nonlin,global_multi_tests,tr_global_multi_tests,global_multi_stats_lin)
%Private, auxiliar file, displaying a global statistical report for multiple
%discriminant analyses. Display will depend wether or not the "linear"
%results are provided for comparison
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%#@ebb: 19,26 Mean (stderr) Jensen-Shannon divergence is displayed and returned
%
%
%__________________________________________________________________________
if nargin>4
    disp('    #Original space of neural activity:'),
    disp(['       Mean missclasified (across validation sets)[SEM]=',num2str(global_multi_stats_lin(1)),'[',num2str(global_multi_stats_lin(2)),']%']),
    %disp(['         [Mean outliers (across validation sets)[SEM]=',num2str(global_multi_stats_lin(3)),'[',num2str(global_multi_stats_lin(4)),']%']),
    disp(['       Mean divergent trajec.(across validation sets)[SEM]=',num2str(global_multi_stats_lin(5)),'[',num2str(global_multi_stats_lin(6)),']%']),
    %disp(['         [Mean outliner trajec. (across validation sets)[SEM]=',num2str(global_multi_stats_lin(7)),'[',num2str(global_multi_stats_lin(8)),']%']),
    disp(['       Mean J-S divergence (across validation sets)[SEM]=',num2str(global_multi_stats_lin(9)),'[',num2str(global_multi_stats_lin(19)),']%']),
end
disp(['    #Expanded space of interactions order ',num2str(order),':']),
disp(['       Mean missclasified (across validation sets)[SEM]=',num2str(global_multi_stats_nonlin(1)),'[',num2str(global_multi_stats_nonlin(2)),']%']),
%disp(['         [Mean outliers (across validation sets)[SEM]=',num2str(global_multi_stats_nonlin(3)),'[',num2str(global_multi_stats_nonlin(4)),']%']),
disp(['       Mean divergent trajec. (across validation sets)[SEM]=',num2str(global_multi_stats_nonlin(5)),'[',num2str(global_multi_stats_nonlin(6)),']%']),
%disp(['         [Mean outliner trajec. (across validation sets)[SEM]=',num2str(global_multi_stats_nonlin(7)),'[',num2str(global_multi_stats_nonlin(8)),']%']),
disp(['       Mean J-S divergence (across validation sets)[SEM]=',num2str(global_multi_stats_nonlin(9)),'[',num2str(global_multi_stats_nonlin(10)),']%']),
disp(' ')
%
if ~all(all(isnan(global_multi_tests))),
    %Misscla.(%)
    disp('    #Normality tests global missclassified (%). p(df,statistic)'),
    disp('       Chi2            Kolmog.-Smirnov            Lilliefords')
    disp(['      ',num2str(global_multi_tests(1,1)),'(',num2str(global_multi_tests(2,1)),',',num2str(global_multi_tests(3,1)),')',...
        '      ',num2str(global_multi_tests(1,2)),'(',num2str(global_multi_tests(2,2)),',',num2str(global_multi_tests(3,2)),')',...
        '      ',num2str(global_multi_tests(1,3)),'(',num2str(global_multi_tests(2,3)),',',num2str(global_multi_tests(3,3)),')']),
    %Comparisons between O=1 and O=order
    if nargin>4
        disp(['   T-test(O=1,O). p(df,t):',num2str(global_multi_tests(1,4)),'(',num2str(global_multi_tests(2,4)),',',num2str(global_multi_tests(3,4)),')']),
        disp(['   Wilcoxon rank-sum(O=1,O). p(n,ranksum):',num2str(global_multi_tests(1,5)),'(',num2str(global_multi_tests(2,5)),',',num2str(global_multi_tests(3,5)),')']),
    end
    %Div.trajec. (%)
    disp('    #Normality tests global divergent trajec.(%). p(df,statistic):'),
    disp('       Chi2            Kolmog.-Smirnov            Lilliefords')
    disp(['      ',num2str(tr_global_multi_tests(1,1)),'(',num2str(tr_global_multi_tests(2,1)),',',num2str(tr_global_multi_tests(3,1)),')',...
        '      ',num2str(tr_global_multi_tests(1,2)),'(',num2str(tr_global_multi_tests(2,2)),',',num2str(tr_global_multi_tests(3,2)),')',...
        '      ',num2str(tr_global_multi_tests(1,3)),'(',num2str(tr_global_multi_tests(2,3)),',',num2str(tr_global_multi_tests(3,3)),')']),
    %Comparisons between O=1 and O=order
    if nargin>4
        disp(['   T-test(O=1,O). p(df,t):',num2str(tr_global_multi_tests(1,4)),'(',num2str(tr_global_multi_tests(2,4)),',',num2str(tr_global_multi_tests(3,4)),')']),
        disp(['   Wilcoxon rank-sum(O=1,O). p(n,ranksum):',num2str(tr_global_multi_tests(1,5)),'(',num2str(tr_global_multi_tests(2,5)),',',num2str(tr_global_multi_tests(3,5)),')']),
    end
end