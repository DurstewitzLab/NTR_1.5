function display_trial(data,labels,trial,order,index_plot)
%feb 011
%WARNING: Private, auxiliary file. Plotting trajectory of a single trial. 
%Different colours correspond to
%different lables -or different phases of the task-
%Hard-coded parameter limiting the number of groups displayed in trajectories. 
%Can be however easily increased by manipulating the code below
%
z_score=2;
max_groups_trajec_display=8;
m_labels=max(labels);l_labels=length(labels);min_labels=min(labels);
n_labels=m_labels-min_labels+1;
%
if n_labels>max_groups_trajec_display,        
        warning('display_trial:maxTr',['Maximum number of groups displayed will be ',num2str(max_groups_trajec_display),' (hard-coded, easy to be changed)'])
end
%
%Check out that all labels are present, either thrown an error. create
%legend text.
%legend_text={};
for i=min_labels:m_labels
    if all(i.*ones(l_labels,1)-labels)
        error('Labels are not natural number increasing by "1". Please create consecutive labels in previous-to-last column of input data'),
    end
    %legend_text=[legend_text,['Task-phase ',num2str(i)]];
end
%
%Centering and normalizing data for better visualization
[n,m]=size(data);
if ~(n==l_labels), error('Labels to be displayed does not match the number of patterns'), end 
minimum=min(min(data));
maximum=max(max(data));
total_max=abs(minimum)+maximum;
centered=(data+abs(minimum).*ones(n,m))./(total_max);
%
%Outliers removal
m_c=mean(centered); std_c=std(centered); smoothed=[];  smoothed_lab=[];
z_up=m_c+z_score*std_c;  z_down=m_c-z_score*std_c;
for i=1:l_labels
        outliers=find(  (centered(i,:)<z_up)&(centered(i,:)>z_down)  );
        if outliers,  
            smoothed=[smoothed;centered(i,:)]; 
            smoothed_lab=[smoothed_lab;labels(i)];
        end
end
disp(['After normalization: ',num2str(length(centered)-length(smoothed)),' outliners deleted ']),
centered=smoothed; labels=smoothed_lab;
[n,~]=size(centered);
%
i=1; count=1;x_prev=centered(1,1);y_prev=centered(1,2);z_prev=centered(1,3);
%
%Next variables are just for plotting the legend. Please see code below.
low_lim=-1.*ones(1,3); up_lim=ones(1,3); 
for i=1:3,
    low_lim(i)=min(centered(:,i));
    up_lim(i)=max(centered(:,i));
end
step=(up_lim(3)-low_lim(3))/n_labels;
margin_1=0.2*(up_lim(1)-low_lim(1));
up_lim(1)=up_lim(1)-margin_1;%Leaving space for z-axis of the rigth plot
flags=ones(1,n_labels);
%
while count<n,
    if labels(count)==1
        x=[x_prev];y=[y_prev];z=[z_prev];
        while (i<=n)&&(labels(i)==1),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'bo','Line','-','LineWidth',1,'MarkerSize',1.5, 'Color',[0,0,.7]);
        title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
         xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        disp('     Displaying task-phase 1...');  
        if flags(1),
           title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
           xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
           text(up_lim(1),low_lim(2)/2,up_lim(3),'Phase 1','Color',[0,0,.7],'FontName','Arial','FontSize',8);
           disp('     Displaying task-phase 1...');
           flags(1)=0;
       end
    elseif (labels(count)==2),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while (i<=n)&&(labels(i)==2),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'ro','Line','-','LineWidth',1,'MarkerSize',1.5, 'Color',[.7,0,0]);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        if flags(2),
           title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
           xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
           text(up_lim(1),low_lim(2)/2,up_lim(3)-step,'Phase 2','Color',[.7,0,0],'FontName','Arial','FontSize',8);
           disp('     Displaying task-phase 2...');
           flags(2)=0;
       end
    elseif (labels(count)==3)&&(n_labels>2),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while(i<=n)&&(labels(i)==3),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'go','Line','-','LineWidth',1,'MarkerSize',1.5, 'Color',[0,.7,0]);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        if flags(3),
           title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
           xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
           text(up_lim(1),low_lim(2)/2,up_lim(3)-2*step,'Phase 3','Color',[0,.7,0],'FontName','Arial','FontSize',8);
           disp('     Displaying task-phase 3...');
           flags(3)=0;
       end
    elseif (labels(count)==4)&&(n_labels>3),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while(i<=n)&&(labels(i)==4)
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'ko','Line','-','LineWidth',1,'MarkerSize',1.5, 'Color',[.4,.4,.4]);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        if flags(4),
           title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
           xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
           text(up_lim(1),low_lim(2)/2,up_lim(3)-3*step,'Phase 4','Color',[.4,.4,.4],'FontName','Arial','FontSize',8);
           disp('     Displaying task-phase 4...');
           flags(4)=0;
       end
    elseif (labels(count)==5)&&(n_labels>4),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while(i<=n)&&(labels(i)==5)
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'yo','Line','-','LineWidth',1,'MarkerSize',1.5,'Color',[.7,.7,0]);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        if flags(5),
           title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
           xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
           text(up_lim(1),low_lim(2)/2,up_lim(3)-4*step,'Phase 4','Color',[.7,.7,0],'FontName','Arial','FontSize',8);
           disp('     Displaying task-phase 5...');
           flags(5)=0;
       end
    elseif (labels(count)==6)&&(n_labels>5),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while(i<=n)&&(labels(i)==6)
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'co','Line','-','LineWidth',1,'MarkerSize',1.5);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        if flags(6),
           title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]),
           xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
           text(up_lim(1),low_lim(2)/2,up_lim(3)-5*step,'Phase 6','Color','c','FontName','Arial','FontSize',8);
           disp('     Displaying task-phase 6...');
           flags(6)=0;
       end
    elseif (labels(count)==7)&&(n_labels>6),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while(i<=n)&&(labels(i)==7),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'mo','Line','-','LineWidth',1,'MarkerSize',1.5);%mo
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        disp('     Displaying task-phase 7...');
    elseif (labels(count)==8)&&(n_labels>7),
        x=[x_prev];y=[y_prev];z=[z_prev];
        while(i<=n)&&(labels(i)==8),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        plot3(x,y,z,'Color',[0.1,0.1,0.1],'Line','-','LineWidth',1,'Marker','o',...
            'MarkerSize',1.5);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        disp('     Displaying task-phase 8...');
    elseif labels(count)==-1
        x=[x_prev];y=[y_prev];z=[z_prev];
        while (i<=n)&&(labels(i)==-1),
            i=i+1;
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
        end
        %
        %Plot in ligth gray, but not during delay period
        if all(labels(count:i,end-1)~=2),
         subplot(2,2,index_plot)  
         plot3(x,y,z,'Color',[0.8,0.8,0.8],'Line','--','Marker','None',...
            'MarkerSize',1);     
        hold on
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        %
    else
        error(['Task-phase number not recognized. It has to be either -1 (no any phase) or ',...
            'a natural number from 1 to ',num2str(m_labels)]),
    end
end
% leg=legend(legend_text);
% set(leg,'FontName','Arial','FontSize',8);
% set(gca,'FontName','Arial','FontSize',8);
xlim([low_lim(1) up_lim(1)]);ylim([low_lim(2) up_lim(2)]);zlim([low_lim(3) up_lim(3)]);
set(gca,'FontName','Arial','FontSize',8);
hold off

