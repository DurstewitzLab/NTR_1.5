function[global_multiv_CE,global_multiv_DIV,global_multiv_JS,global_multiv_L,global_multiv_L2,global_multiv_M,...
    twoClass_CE,twoClass_DIV,twoClass_JS,twoClassL,twoCLassM,twoClassL2,...
     multiv_classCE, multiv_classDIV, multiv_classJS]=bootstrap_generation(file_name,conf,n_boot)
%Auxiliar batch function for surrogate generation
laboratories=1;
if nargin<3,
    n_boot=50;
end
if nargin<2,
    conf=setup_config(1,0.148,1,0,0,0,[],1,2,0,15,1,30,1,30);
end
if nargin<1,
    file_name='Saline_1mg_A1_11_12_14_19';
end
try
    matlabpool (laboratories); disp(['Engaging...',num2str(matlabpool ('size')),' labs']),
catch OPENSESSION
    disp('Closing previous sessions...'); matlabpool close; matlabpool (laboratories); disp([num2str(matlabpool ('size')),' labs involved']),
end
%
global_multiv_CE=zeros(n_boot,1); 
global_multiv_DIV=zeros(n_boot,1);
global_multiv_JS=zeros(n_boot,1);
%
multiv_classCE=zeros(n_boot,6);
multiv_classDIV=zeros(n_boot,6);
multiv_classJS=zeros(n_boot,6);
%
twoClass_CE=zeros(n_boot,15);
twoClass_DIV=zeros(n_boot,15);
twoClass_JS=zeros(n_boot,15);%Strange 
%
%twoClassJS=zeros(n_boot,1);
twoClassL=zeros(n_boot,1);
twoCLassM=zeros(n_boot,1);
twoClassL2=zeros(n_boot,1);

%parfor
for bot=1:n_boot
    disp(' ')
    disp('***************************************')
    disp('***'), disp(['Boot. ',num2str(bot),'/',num2str(n_boot)]),disp('***')
    disp('***************************************')
    disp(' ')
    results=ntr(file_name,conf);
    CE_multiv=results.multiClass_err; 
    DE_multiv=results.multiClass_diver;
    JS=results.js;
    L=results.lambda;
    L2=results.lambda2;
    M=results.mhd;
    %
    class_multi=results.cum_multi_stats;
    multiv_classCE(bot,:)=mean(class_multi(:,:,1));
    multiv_classDIV(bot,:)=mean(class_multi(:,:,2));
    multiv_classJS(bot,:)=mean(class_multi(:,:,3));
    %
    pair_CE_two=results.twoClass_mean_err; 
    pair_DE_two=results.twoClass_mean_div;
    pair_JS_two=results.twoClass_mean_js;
    
   %js_2=results.twoClass_2mean_js;
    l_2=results.twoClass_2mean_lambda;
    m_2=results.twoClass_2mean_mhd;
    l2_2=results.twoClass_2mean_lambda2; 
    %
    global_multiv_CE(bot)=CE_multiv;
    global_multiv_DIV(bot)=DE_multiv;
    global_multiv_JS(bot)=JS;
    global_multiv_L(bot)=L;
    global_multiv_L2(bot)=L2;
    global_multiv_M(bot)=M;
    %
    twoClass_CE(bot,:)=pair_CE_two';
    twoClass_DIV(bot,:)=pair_DE_two';
    twoClass_JS(bot,:)=pair_JS_two';
    twoClassL(bot)= l_2;
    twoCLassM(bot)=m_2;
    twoClassL2(bot)=l2_2;  
    
end
matlabpool close;