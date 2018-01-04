function generateLabels(data1,e1,data2,e2)
n_epochs=6;
min_epoch=1;
combi=combnk(1:n_epochs,2); n_combi=length(combi(:,1));
combi=combi+(min_epoch-1).*ones(n_combi,2);     
XTickLabels=cell(1,n_combi);
for i=1:n_combi,
        current_pair=combi(i,:);
        if strcmp(num2str(current_pair(1)),'1')
            s1='TrC';
        elseif  strcmp(num2str(current_pair(1)),'2')
            s1='TrR';
        elseif  strcmp(num2str(current_pair(1)),'3')
            s1='D';
        elseif  strcmp(num2str(current_pair(1)),'4')
            s1='TsC';
        elseif  strcmp(num2str(current_pair(1)),'5')
            s1='TsR';
        elseif  strcmp(num2str(current_pair(1)),'6')
            s1='I';
        end
        if strcmp(num2str(current_pair(2)),'1')
            s2='TrC';
        elseif  strcmp(num2str(current_pair(2)),'2')
            s2='TrR';
        elseif  strcmp(num2str(current_pair(2)),'3')
            s2='D';
        elseif  strcmp(num2str(current_pair(2)),'4')
            s2='TsC';
        elseif  strcmp(num2str(current_pair(2)),'5')
            s2='TsR';
        elseif  strcmp(num2str(current_pair(2)),'6')
            s2='I';
        end    
        XTickLabels{i}=[s1,'-',s2];
end
display_class_results(data1,e1,data2,e2,4,'',XTickLabels)
%
index_Trc=[1:5];index_Trr=[1,6,7,8,9];index_D=[2,6,10,11,12];
index_Tsc=[3,7,10,13,14];index_Tsr=[4,8,11,13,15];index_I=[5,9,12,14,15];
%
Trc=mean(data1(index_Trc));e_Trc=mean(e1(index_Trc));
Trr=mean(data1(index_Trr));e_Trr=mean(e1(index_Trr));
D=mean(data1(index_D));e_D=mean(e1(index_D));
Tsc=mean(data1(index_Tsc));e_Tsc=mean(e1(index_Tsc));
Tsr=mean(data1(index_Tsr));e_Tsr=mean(e1(index_Tsr));
I=mean(data1(index_I));e_I=mean(e1(index_I));
%
mdata1=[Trc,Trr,D,Tsc,Tsr,I]'; 
e_m1=[e_Trc,e_Trr,e_D,e_Tsc,e_Tsr,e_I]'; 
%
Trc=mean(data2(index_Trc));e_Trc=mean(e2(index_Trc));
Trr=mean(data2(index_Trr));e_Trr=mean(e2(index_Trr));
D=mean(data2(index_D));e_D=mean(e2(index_D));
Tsc=mean(data2(index_Tsc));e_Tsc=mean(e2(index_Tsc));
Tsr=mean(data2(index_Tsr));e_Tsr=mean(e2(index_Tsr));
I=mean(data2(index_I));e_I=mean(e2(index_I));
%
mdata2=[Trc,Trr,D,Tsc,Tsr,I]'; 
e_m2=[e_Trc,e_Trr,e_D,e_Tsc,e_Tsr,e_I]'; 
%
disp('****Comparisons ttests***')
[h,p,ci,s]=ttest2(data1(index_Trc),data2(index_Trc));disp(['TrC(p)=',num2str(p)]),
s
[h,p,ci,s]=ttest2(data1(index_Trr),data2(index_Trr));disp(['Trr(p)=',num2str(p)]),
s
[h,p,ci,s]=ttest2(data1(index_D),data2(index_D));disp(['D(p)=',num2str(p)]),
s
% [h,p] = lillietest(data2(index_D))
% [h,p] = lillietest(data1(index_D))
% [h,p,ksstat,cv] = kstest(data1(index_D))
% [h,p,ksstat,cv] = kstest(data2(index_D))
[h,p,ci,s]=ttest2(data1(index_Tsc),data2(index_Tsc));disp(['Tsc(p)=',num2str(p)]),
s
% [h,p] = lillietest(data2(index_Tsc))
% [h,p] = lillietest(data1(index_Tsc))
[h,p,ci,s]=ttest2(data1(index_Tsr),data2(index_Tsr));disp(['Tsr(p)=',num2str(p)]),

s
[h,p,ci,s]=ttest2(data1(index_I),data2(index_I));disp(['I(p)=',num2str(p)]),
s
disp('***********************')
disp('****Comparisons ranskum***')
[p,h,s]=ranksum(data1(index_Trc),data2(index_Trc));disp(['TrC(p)=',num2str(p)]),
s
[p,h,s]=ranksum(data1(index_Trr),data2(index_Trr));disp(['Trr(p)=',num2str(p)]),
s
[p,h,s]=ranksum(data1(index_D),data2(index_D));disp(['D(p)=',num2str(p)]),
s
[p,h,s]=ranksum(data1(index_Tsc),data2(index_Tsc));disp(['Tsc(p)=',num2str(p)]),
s
[p,h,s]=ranksum(data1(index_Tsr),data2(index_Tsr));disp(['Tsr(p)=',num2str(p)]),
s
[p,h,s]=ranksum(data1(index_I),data2(index_I));disp(['I(p)=',num2str(p)]),
s
display_class_results(mdata1,e_m1,mdata2,e_m2,4,'averaged');
