%Batch procces script
%
%I-Config
clear; clc; close all;
%Windows
%save_dir='../ntr1_5a/results/';
%Linux
save_dir='/home/emili.balaguer/ntr1_5b/results/';
n_boot=1;
%n_boot=2;
%params=orders in this case
param_name='order';
 %array_vals=[1,2,3,4,5,6,8,10];%One parameter will be changed in turn and bootstrap
%                        will be generated for each
array_vals=[1];
%array_vals=[1,3,5];
n=length(array_vals);
%file_name='Saline_1mg_A1_11_12_14_19';
%**Ah-hoc-setting***
tau=.1;
%tau=.2;
data_lim=30;
%data_lim=5;
file_name='Amph_1mg_A1_11_12_14_19';
% file_name='Saline_3mg_A8_9_13_16_20.mat';
%file_name='Amph_3mg_A8_9_13_16_20.mat';
%file_name='Saline_3mg_8_9_13_16_17_20.mat';
%file_name='Amph_3mg_8_9_13_16_17_20.mat';
%
%"Basic" configuration
%basic_conf=setup_config(1,0.148,0,0,0,0,[],1,2,0,13,1,30,1,30);%Basic configuration for saline 1mg
%basic_conf=setup_config(1,0.41,0,0,0,0,[],1,2,0,10,1,30,1,data_lim);%Basic configuration for amph 1mg
%conf=setup_config(5,0.148,1,0,0,0,[],1,2,0,15,1,30,1,30);
basic_conf=setup_config(1,0.423,0,0,0,0,[],1,2,0,13,1,30,1,data_lim);
%
%We will open a whole set of laboratories for each
%
%II-Process (means and standard errors).
%First two columns contain non-bootstraped data
params_multi_CE=zeros(n,4);
params_multi_DIV=zeros(n,4);
params_multi_JS=zeros(n,4);
params_multi_lambda=zeros(n,4);
params_multi_lambda2=zeros(n,4);
params_multi_mhd=zeros(n,4);
%
params_multi_class_CE=[];%Means per class
params_multi_class_DIV=[];
params_multi_class_JS=[];
%Boot
Bparams_multi_class_CE=[];%Means per class
Bparams_multi_class_DIV=[];
Bparams_multi_class_JS=[];
%
params_2class_CE=zeros(n,4);
params_2class_DIV=zeros(n,4);
params_2class_JS=zeros(n,4);
params_2class_lambda=zeros(n,4);
params_2class_lambda2=zeros(n,4);
params_2class_mhd=zeros(n,4);
%
params_2_CE=[];
params_2_DIV=[];
params_2_JS=[];
Bparams_2_CE=[];
Bparams_2_DIV=[];
Bparams_2_JS=[];
%
for i=1:n,
    disp(''),
    disp('_______________________'),
    disp('_______________________'),
    disp(''),
    disp([num2str(array_vals(i)),'-th ',param_name]),
    disp(''),
    tic
    %
    %Setting the specific configuration of this parameter.
    %Now, the paramter is just the order (can be easily modfied)
    conf=basic_conf;
    conf.order=array_vals(i);
    %     % STARTS COMMENTS
    %     %II.1 Non-bootstraped data
    results=ntr(file_name,conf);
    params_multi_CE(i,1)=results.multiClass_err;
    params_multi_CE(i,2)= results.bars_multiClass_err;
    params_multi_DIV(i,1)=results.multiClass_diver;
    params_multi_DIV(i,2)= results.bars_multiClass_diver;
    params_multi_JS(i,1)=results.js;
    params_multi_JS(i,2)=results.bars_js;
    params_multi_lambda(i,1)=results.lambda;
    params_multi_lambda(i,2)=results.bars_lambda;
    params_multi_lambda2(i,1)=results.lambda2;
    params_multi_lambda2(i,2)=results.bars_lambda2;
    params_multi_mhd(i,1)=results.mhd;
    params_multi_mhd(i,2)=results.bars_mhd;
    %
    disp(['Original data: Mean multiv.CE(STDERR): ',num2str(params_multi_CE(i,1)),'(',num2str(params_multi_CE(i,2)),')']);
    disp(['Original data: Mean multiv. DIV(STDERR): ',num2str(params_multi_DIV(i,1)),'(',num2str(params_multi_DIV(i,2)),')']);
    %
    class_multi=results.cum_multi_stats;
    if n<2
        params_multi_class_CE=[mean(class_multi(:,:,1))',std(class_multi(:,:,1))'];%Means and errors per class
        params_multi_class_DIV=[mean(class_multi(:,:,2))',std(class_multi(:,:,2))'];
        params_multi_class_JS=[mean(class_multi(:,:,3))',std(class_multi(:,:,3))'];
    else
        warning('meta_boot_batch:no_class_results','Warning class resutls only supported for a single order')
    end
    %
    %2-class. Here is different
    params_2class_CE(i,1)=results.twoClass_2mean_err;
    params_2class_CE(i,2)= results.twoClass_2std_err;
    params_2class_DIV(i,1)=results.twoClass_2mean_div;
    params_2class_DIV(i,2)=results.twoClass_2std_div;
    params_2class_JS(i,1)=results.twoClass_2mean_js;
    params_2class_JS(i,2)=results.twoClass_2std_js;
    params_2class_lambda(i,1)=results.twoClass_2mean_lambda;
    params_2class_lambda(i,2)=results.twoClass_2std_lambda;
    params_2class_lambda2(i,1)=results.twoClass_2mean_lambda2;
    params_2class_lambda2(i,2)=results.twoClass_2std_lambda2;
    params_2class_mhd(i,1)=results.twoClass_2mean_mhd;
    params_2class_mhd(i,2)=results.twoClass_2std_mhd;
    
    %Recoverign all-combi-results
    
    CE2=results.twoClass_mean_err;
    eCE2=results.twoClass_std_err;
    DIV2=results.twoClass_mean_div;
    eDIV2=results.twoClass_std_div;
    JS2=results.twoClass_mean_js;
    eJS2=results.twoClass_std_js;
    
    disp(['Original data: Mean 2-class (STDERR): ',num2str(params_2class_CE(i,1)),'(',num2str(params_2class_CE(i,2)),')']);
    disp(['Original data: Mean 2-class DIV(STDERR): ',num2str(params_2class_DIV(i,1)),'(',num2str(params_2class_DIV(i,2)),')']);
    
    if n<2
        params_2_CE=[CE2,eCE2];%Means and errors per class
        params_2_DIV=[DIV2,eDIV2];
        params_2_JS=[JS2,eJS2];
    else
        warning('meta_boot_batch:no_class_results','Warning class results only supported for a single order')
    end
    %
    %II.2 Bootstraped data
    conf.shuff_epoch=1;
    [global_multiv_CE,global_multiv_DIV,global_multiv_JS,global_multiv_L,global_multiv_L2,global_multiv_M,...
        twoClass_CE,twoClass_DIV,twoClass_JS,twoClassL,twoClassM,twoClassL2,...
        multiv_classCE, multiv_classDIV, multiv_classJS]=bootstrap_generation(file_name,conf,n_boot);
    %
    %Multivariate
    mean_multi_CE=mean(global_multiv_CE);
    err_multi_CE=std(global_multiv_CE)/sqrt(length(global_multiv_CE));
    mean_multi_DIV=mean(global_multiv_DIV);
    err_multi_DIV=std(global_multiv_DIV)/sqrt(length(global_multiv_DIV));
    mean_multi_JS=mean(global_multiv_JS);
    err_multi_JS=std(global_multiv_JS)/sqrt(length(global_multiv_JS));
    mean_multi_lambda=mean(global_multiv_L);
    err_multi_lambda=std(global_multiv_L)/sqrt(length(global_multiv_L));
    mean_multi_lambda2=mean(global_multiv_L2);
    err_multi_lambda2=std(global_multiv_L2)/sqrt(length(global_multiv_L2));
    mean_multi_mhd=mean(global_multiv_M);
    err_multi_mhd=std(global_multiv_M)/sqrt(length(global_multiv_M));
    %
    disp([num2str(n_boot),' Bootstraps of file ',file_name,':']),
    disp(['Mean multiv.CE(STDERR): ',num2str(mean_multi_CE),'(',num2str(err_multi_CE),')']);
    disp(['Mean multiv. DIV(STDERR): ',num2str(mean_multi_DIV),'(',num2str(err_multi_DIV),')']);
    
    %
    %Pairwise
    %Averaging across bootstraps and considering this averaged series
    twoClass_CE=mean( twoClass_CE); twoClass_DIV=mean( twoClass_DIV);twoClass_JS=mean( twoClass_JS);
    %
    mean_two_CE=mean(twoClass_CE);
    err_two_CE=std(twoClass_CE)/sqrt(length(twoClass_CE));
    mean_two_DIV=mean(twoClass_DIV);
    err_two_DIV=std(twoClass_DIV)/sqrt(length(twoClass_DIV));
    mean_two_JS=mean(twoClass_JS);
    err_two_JS=std(twoClass_JS)/sqrt(length(twoClass_JS));
    %
    disp(['Mean 2-class CE(STDERR): ',num2str(mean_two_CE),'(',num2str(err_two_CE),')']);
    disp(['Mean 2-class DIV(STDERR): ',num2str(mean_two_DIV),'(',num2str(err_two_DIV),')']);
    %
    %         mean_two_JS=mean( twoClassJS);
    %         err_two_JS=std(twoClassJS)/sqrt(length(twoClassJS));
    mean_two_L=mean( twoClassL);
    err_two_L=std(twoClassL)/sqrt(length(twoClassL));
    mean_two_M=mean( twoClassM);
    err_two_M=std(twoClassM)/sqrt(length(twoClassM));
    mean_two_L2=mean( twoClassL2);
    err_two_L2=std(twoClassL2)/sqrt(length(twoClassL2));
    %
    %
    h=floor(toc/3600); m=floor(toc/60)-(h*60);s=toc-(h*3600)-(m*60);
    disp(['Total time ',num2str(h),':',num2str(m),':',num2str(s)]),
    %
    params_multi_CE(i,3)=mean_multi_CE;
    params_multi_CE(i,4)= err_multi_CE;
    params_multi_DIV(i,3)=mean_multi_DIV;
    params_multi_DIV(i,4)= err_multi_DIV;
    params_multi_JS(i,3)=mean_multi_JS;
    params_multi_JS(i,4)=err_multi_JS;
    params_multi_lambda(i,3)=mean_multi_lambda;
    params_multi_lambda(i,4)=err_multi_lambda;
    params_multi_lambda2(i,3)=mean_multi_lambda2;
    params_multi_lambda2(i,4)=err_multi_lambda2;
    params_multi_mhd(i,3)=mean_multi_mhd;
    params_multi_mhd(i,4)=err_multi_mhd;
    %
    if n<2
        Bparams_multi_class_CE=multiv_classCE;
        Bparams_multi_class_DIV=multiv_classDIV;
        Bparams_multi_class_JS=multiv_classJS;
    end
    %
    params_2class_CE(i,3)=mean_two_CE;
    params_2class_CE(i,4)= err_two_CE;
    params_2class_DIV(i,3)=mean_two_DIV;
    params_2class_DIV(i,4)= err_two_DIV;
    params_2class_JS(i,3)=mean_two_JS;
    params_2class_JS(i,4)=err_two_JS;
    params_2class_lambda(i,3)= mean_two_L;
    params_2class_lambda(i,3)=err_two_L;
    params_2class_lambda2(i,3)=mean_two_L2;
    params_2class_lambda2(i,4)=err_two_L2;
    params_2class_mhd(i,3)=mean_two_M;
    params_2class_mhd(i,4)= err_two_M;
    %
    if n<2
        Bparams_2_CE=twoClass_CE;
        Bparams_2_DIV=twoClass_DIV;
        Bparams_2_JS=twoClass_JS;
    end
    %
end
%III. Saving and displaying
if n>1
name=['Bootstraps_',file_name,'_','Regul',num2str(conf.regul),'_Tau',num2str(tau),'_n',...
    num2str(n_boot),'_DataLim',num2str(data_lim),'.mat'];
else
    name=['Bootstraps_',file_name,'_Order',num2str(array_vals(1)),'_Regul',num2str(conf.regul),'_Tau',num2str(tau),'_n',...
    num2str(n_boot),'_DataLim',num2str(data_lim),'.mat'];
end
save([save_dir,name]);
disp(['Results saved in: ',which(name)]),
warning('on','all')
display_orders_bootstraps(params_multi_CE,params_multi_DIV,params_2class_CE,params_2class_DIV,array_vals,name);


