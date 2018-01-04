function basic_stats=fill_predic(basic_stats,set_lab_ref,set_lab)
%Private, auxiliar file for cleaning code. Do not modify.
n_epochs=length(set_lab_ref);
%Real values
true_se=basic_stats.cum_misscla_class;
true_bins=basic_stats.n_bins_per_class;
true_div=basic_stats.cum_diverg_trajec_class;
true_tr=basic_stats.cum_trajec_class;
true_prob=basic_stats.posterior_probs_vectors;
true_lik=basic_stats.lik_vectors;
true_class_js=basic_stats.class_js;
n_patterns=length(true_prob(1,:));
%
%Next vectors will be filled with real or averaged predictions
fill_se_class=zeros(1,n_epochs);fill_bins_class=zeros(1,n_epochs);
fill_diverg_trajec_class=zeros(1,n_epochs);fill_trajec_class=zeros(1,n_epochs);
fill_prob=zeros(n_epochs,n_patterns); fill_lik=zeros(n_epochs,n_patterns);
fill_class_js=zeros(1,n_epochs);
%
count=1;
for i=1:length(set_lab_ref),
    if any(set_lab==set_lab_ref(i)),
        fill_se_class(i)=true_se(count);
        fill_bins_class(i)=true_bins(count);
        fill_diverg_trajec_class(i)=true_div(count);
        fill_trajec_class(i)=true_tr(count);
        fill_prob(i,:)=true_prob(count,:);
        fill_lik(i,:)=true_lik(count,:);
         fill_class_js(i)=true_class_js(count);
        count=count+1;
    else
        fill_se_class(i)=mean(true_se);
        fill_bins_class(i)=mean(true_bins);
        fill_diverg_trajec_class(i)=mean(true_div);
        fill_trajec_class(i)=mean(true_tr);
        fill_prob(i,:)=zeros(1,n_patterns);
        fill_lik(i,:)=zeros(1,n_patterns);
        fill_class_js(i)=min(true_class_js);
    end
end
basic_stats.cum_misscla_class=fill_se_class;
basic_stats.n_bins_per_class=fill_bins_class;
basic_stats.cum_diverg_trajec_class=fill_diverg_trajec_class;
basic_stats.cum_trajec_class=fill_trajec_class;
basic_stats.posterior_probs_vectors=fill_prob;
basic_stats.lik_vectors=fill_lik;
basic_stats.class_js=fill_class_js;