function test_results=stat_tests(data_nonlin,data_lin)
%auxiliar file performing a battery of Statistical comparisons.
%For a single input variable, only nomality tests will be performed.
%
% For a single input
%   Column 1: Row 1: p-value of Chi2gof. Row 2: df. Row 3: stastitic value.
%   Column 2: Row 1: p-value of ks. Row 2: n. Row 3: stastitic value.
%   Column 3: Row 1: p-value of Liliefors. Row 2: n. Row 3: stastitic value.
%
% For two inputs
% "test_results" matrix will have the next content:
%   Column 1: Row 1: mimimum p-value of Chi2gof(O=1,O). Row 2: df. Row 3: stastitic value.
%   Column 2: mimimum p-value of ks(O=1,O). Row 2: n. Row 3: stastitic value.
%   Column 3: mimimum p-value of Liliefors(O=1,O). Row 2: n. Row 3: stastitic value.
%   Column 4: p-value of ttest(O=1,O). Row 2: df. Row 3: stastitic value.
%   Column 5: p-value of ranksum(O=1,O). Row 2: df. Row 3: stastitic value
if nargin<2
    Nans=isnan(data_nonlin);
    if length(Nans(Nans>0))<length(data_nonlin)/2,
        data_nonlin=data_nonlin(~isnan(data_nonlin));
        [~,p,stats]=chi2gof(data_nonlin);
        test_results(1,1)=p;test_results(2,1)=stats.df;test_results(3,1)=stats.chi2stat;
        %
        %Note: distribution has to be standarized for the kolmog-smir. test
        [~,p,ksstat]=kstest(zscore(data_nonlin));
        test_results(1,2)=p;test_results(2,2)=length(data_nonlin);test_results(3,2)=ksstat;
        %
        [~,p,lilstat]=lillietest(zscore(data_nonlin));
        test_results(1,3)=p;test_results(2,3)=length(data_nonlin);test_results(3,3)=lilstat;
    else
        test_results=NaN.*ones(3,3);
    end
else
    Nans1=isnan(data_nonlin);Nans2=isnan(data_lin);
    if (length(Nans1(Nans1>0))<length(data_nonlin)/2)&&(length(Nans2(Nans2>0))<length(data_lin)/2),
        data_lin=data_lin(~isnan(data_lin));  data_nonlin=data_nonlin(~isnan(data_nonlin));
        %a) Normality tests (p<0.05 implies departure from gaussian)
        [~,p1,stats1]=chi2gof(data_lin);[~,p2,stats2]=chi2gof(data_nonlin);
        if p1>p2, p=p2; df=stats2.df; stat=stats2.chi2stat; else p=p1; df=stats1.df; stat=stats1.chi2stat; end
        test_results(1,1)=p;test_results(2,1)=df;test_results(3,1)=stat;
        %
        %Note: distribution has to be standarized for the kolmog-smir. test
        
        [h1,p1,stats1]=kstest(zscore(data_lin));[h2,p2,stats2]=kstest(zscore(data_nonlin));
        if p1>p2,
            p=p2; stat=stats2;
        else
            p=p1; stat=stats1;
        end
        test_results(1,2)=p;test_results(2,2)=length(data_lin);test_results(3,2)=stat;
        %
        [~,p1,stats1]=lillietest(zscore(data_lin));[~,p2,stats2]=lillietest(zscore(data_nonlin));
        if p1>p2, p=p2; stat=stats2; else p=p1; stat=stats1; end
        test_results(1,3)=p;test_results(2,3)=length(data_lin);test_results(3,3)=stat;
        %
        %b) Comparisons
        %t-tests (p<0.05 rejects equal means)
        
        [~,p,~,stats]=ttest2(data_lin,data_nonlin,[],[],'unequal');
        test_results(1,4)=p;test_results(2,4)=stats.df;test_results(3,4)=stats.tstat;
        %
        %Two-sided Wilcoxon rank-sum nonparametric test (p<0.05 rejects equal means
        [p,~,stats]=ranksum(data_lin,data_nonlin);
        test_results(1,5)=p;test_results(2,5)=length(data_lin);test_results(3,5)=stats.ranksum;
    else
        test_results=NaN.*ones(3,5);
    end
end
