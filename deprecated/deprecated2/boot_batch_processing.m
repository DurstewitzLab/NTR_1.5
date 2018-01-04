%Batch procces script
%
%I-Config
clear; clc; close all;
%Windows
%save_dir='../ntr1_5a/results/';
%Linux
save_dir='..\results\';
n_boot=100;
laboratories=8;
%file_name='Saline_1mg_A1_11_12_14_19';
% file_name='Amph_1mg_A1_11_12_14_19';
% file_name='Saline_3mg_A8_9_13_16_20.mat';
% file_name='Amph_3mg_A8_9_13_16_20.mat';
% file_name='Saline_3mg_8_9_16_17_20.mat';
 file_name='Amph_3mg_8_9_13_16_17_20.mat';
%
%II-Process
%conf=setup_config(5,0.148,1,0,0,0,[],1,2,0,30,1,30,1,30);
%conf=setup_config(5,0.148,1,0,0,0,[],1,2,0,15,1,30,1,30);
conf=setup_config(1,0.385,1,0,0,0,[],1,2,0,15,1,30,1,30);
try
    matlabpool (laboratories); disp([num2str(matlabpool ('size')),' labs involved']),
catch OPENSESSION
    disp('Closing previous sessions...'); matlabpool close; matlabpool (laboratories); disp([num2str(matlabpool ('size')),' labs involved']),
end
tic
[global_multiv_CE,global_multiv_DIV,twoClass_CE,twoClass_DIV]=bootstrap_generation(file_name,conf,n_boot);
%
mean_multi_CE=mean(global_multiv_CE);
err_multi_CE=std(global_multiv_CE)/sqrt(length(global_multiv_CE));
mean_multi_DIV=mean(global_multiv_DIV);
err_multi_DIV=std(global_multiv_DIV)/sqrt(length(global_multiv_DIV));
disp([num2str(n_boot),' bootstraps of file ',file_name,':']),
disp(['Mean CE(STDERR): ',num2str(mean_multi_CE),'(',num2str(err_multi_CE),')']);
disp(['Mean DIV(STDERR): ',num2str(mean_multi_DIV),'(',num2str(err_multi_DIV),')']);
%
name=['O',num2str(conf.order),'_Regul',num2str(conf.regul),'_Boot',num2str(n_boot),'_',file_name,'.mat'];
save([save_dir,name],'global_multiv_CE','global_multiv_DIV','twoClass_CE','twoClass_DIV');
disp(['Results saved in: ',which(name)]),
matlabpool close;
h=floor(toc/3600); m=floor(toc/60)-(h*60);s=toc-(h*3600)-(m*60);
disp(['Total time ',num2str(h),':',num2str(m),':',num2str(s)]),

