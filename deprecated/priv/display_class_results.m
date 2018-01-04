function display_class_results(data1,e1,data2,e2,order,fig_title,combi)
%Private, auxiliar file, displaying means and SEM per class or
%pairwise comparisons.
%Inputs
%   data1 (n_comparisons x 2) = Each column is the % missclassified vectors
%                              for two specified "orders".
%                              If only one column specified, it correspond
%                              to O="order".
%   e1 (n_comparisons x 2) = Error bars corresponding to data1 columns
%   data2 (n_comparisons x 2) = Same for the % of divergent trajectories
%   e2 (n_comparisons x 2) = Error bars corresponding to data2 columns
%   order (1x1) = Nonlinear expansion order
%   title (string)= Text describing the figure
%

[n_comp,n_orders]=size(data1);%Number of comparisons
XTickLabels=cell(1,n_comp);
for i=1:n_comp,
    if (nargin<7),
        XTickLabels{i}=num2str(i);
    else
        if ~any(isletter(combi))
            XTickLabels=combi;
        else
            current_pair=combi(i,:);
            XTickLabels{i}=[num2str(current_pair(1)),'-',num2str(current_pair(2))];
        end
    end
end
if (length(data2(:,1))~=n_comp)||(length(e1(:,1))~=n_comp)||(length(e2(:,1))~=n_comp),
    warning('display_class_results:data_mismatch','Mismatch between data length, figure not produced'),
else
    figure
    if n_orders>1,
        long_title=[fig_title,': original and expanded space (order=',num2str(order),')'];
    else
        long_title=[fig_title,': expanded space (order=',num2str(order),')'];
    end
        set(gcf,'Color',[1,1,1],'Name',long_title),
        %
        %Missclasified vectors
        subplot(1,2,1),errorbar(data1(:,1),e1(:,1),'bo');
        if n_orders>1, hold on,errorbar(data1(:,2),e2(:,2),'ro'); end
        xl=xlabel('Task-epoch (or task-epoch-pair comparison)');
        yl=ylabel('Segregation error (%)');
        set(gca,'XTick',(1:n_comp),'XTickLabelMode','manual','XTickLabel',XTickLabels,...
            'FontSize',8);
        tit=title([fig_title,': Missclasified (%)']);
        if n_orders>1, 
            leg=legend('Activity space',['Expanded space of order ',num2str(order)]);
            set(leg,'FontName','Arial','FontSize',10);
        end
        set(tit,'FontName','Arial','FontSize',10,'FontAngle','Italic');
        set(xl,'FontName','Arial','FontSize',10);
        set(yl,'FontName','Arial','FontSize',10);
        set(gca,'FontName','Arial','FontSize',8);
        view(90,90);
        %
        %Divergent trajectories
        subplot(1,2,2),errorbar(data2(:,1),e2(:,1),'go');
        if n_orders>1, hold on, errorbar(data2(:,2),e2(:,2),'co'); end
        xl=xlabel('Task-epoch (or task-epoch-pair comparison)');
        yl=ylabel('Divergent trajectories (%)');
        set(gca,'XTick',(1:n_comp),'XTickLabelMode','manual','XTickLabel',XTickLabels,...
            'FontSize',8);
        tit=title([fig_title,': Divergent trajectories (%)']);
        if n_orders>1, 
            leg=legend('Activity space',['Expanded space of order ',num2str(order)]); 
            set(leg,'FontName','Arial','FontSize',10);
        end
        set(tit,'FontName','Arial','FontSize',10,'FontAngle','Italic');
        set(xl,'FontName','Arial','FontSize',10);
        set(yl,'FontName','Arial','FontSize',10);
        set(gca,'FontName','Arial','FontSize',8);
        view(90,90);
end