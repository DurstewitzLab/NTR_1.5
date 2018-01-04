
function class_multi_stats_report(order,cum_multi_stats_nonlin,multi_tests,tr_multi_tests,full_stats_disp,cum_multi_stats_lin)
%Private, auxiliar file, displaying a class-by-class statistical report for multiple
%discriminant analyses. Display will depend wether or not the "linear"
%results are provided for comparison
%
%Figure displaying % misccla. and % of divergent trajec
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
%
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%
disp('    '),
disp('    Displaying % missclassified, % divergent trajec per task-epoch...')
if nargin<4, full_stats_disp=1; end
if nargin>5
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
if  (full_stats_disp>1)
    if (~all(all(isnan(multi_tests))))
        n_epochs=length(cum_multi_stats_nonlin(1,:,1));
        %Misscla.(%).
        disp('    #Iask-epochs misscla. (%). Normality tests. p(df,statistic):'),
        disp('Task-epoch/Chi2          Kolmog.-Smirnov            Lilliefords')
        for i=1:n_epochs
            disp([num2str(i),'      ',num2str(multi_tests(i,1,1)),'(',num2str(multi_tests(i,2,1)),',',num2str(multi_tests(i,3,1)),')',...
                '      ',num2str(multi_tests(i,1,2)),'(',num2str(multi_tests(i,2,2)),',',num2str(multi_tests(i,3,2)),')',...
                '      ',num2str(multi_tests(i,1,3)),'(',num2str(multi_tests(i,2,3)),',',num2str(multi_tests(i,3,3)),')']),
        end
        if nargin>5
            %Class-by class comparisons between O=1 and O=order
            disp('Task-epoch          T-test(O=1,O). p(df,t)          Wilcoxon rank-sum(O=1,O). p(n,ranksum)')
            for i=1:n_epochs
                disp([num2str(i),'       ',num2str(multi_tests(i,1,4)),'(',num2str(multi_tests(i,2,4)),',',num2str(multi_tests(i,3,4)),')',...
                    '       ',num2str(multi_tests(i,1,5)),'(',num2str(multi_tests(i,2,5)),',',num2str(multi_tests(i,3,5)),')']),
            end
        end
        %Div.trajec. (%). Class-by class comparisons between O=1 and O=order
        %Only displayed if "full_stats" is selected.
        disp('     #Task-epochs divergent trajec. (%). Normality tests. p(df,statistic):'),
        disp('Task-epoch/Chi2          Kolmog.-Smirnov            Lilliefords')
        for i=1:n_epochs
            disp([num2str(i),'      ',num2str(tr_multi_tests(i,1,1)),'(',num2str(tr_multi_tests(i,2,1)),',',num2str(tr_multi_tests(i,3,1)),')',...
                '      ',num2str(tr_multi_tests(i,1,2)),'(',num2str(tr_multi_tests(i,2,2)),',',num2str(tr_multi_tests(i,3,2)),')',...
                '      ',num2str(tr_multi_tests(i,1,3)),'(',num2str(tr_multi_tests(i,2,3)),',',num2str(tr_multi_tests(i,3,3)),')']),
        end
        if nargin>5
            disp('Task-epoch          T-test(O=1,O). p(df,t)          Wilcoxon rank-sum(O=1,O). p(n,ranksum)')
            for i=1:n_epochs
                disp([num2str(i),'       ',num2str(tr_multi_tests(i,1,4)),'(',num2str(tr_multi_tests(i,2,4)),',',num2str(tr_multi_tests(i,3,4)),')',...
                    '       ',num2str(tr_multi_tests(i,1,5)),'(',num2str(tr_multi_tests(i,2,5)),',',num2str(tr_multi_tests(i,3,5)),')']),
            end
        end
    end
end


