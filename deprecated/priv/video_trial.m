function video_trial(data,labels,trial,order,index_plot)
%feb 011.
%WARNING: Private, auxiliary file. Used straight for the "nti4_0b" open toolbox,
%Hard-coded parameter limiting the number of groups displayed in trajectories.
%Can be however easily increased by manipulating the code below.
%

max_groups_trajec_display=6;
m_labels=max(labels);l_labels=length(labels);min_labels=min(labels);
n_labels=m_labels-min_labels+1;
%
if n_labels>max_groups_trajec_display,
    warning('video_trial:maxGroups',['Maximum number of groups visualized will be ',num2str(max_groups_trajec_display),' (hard-coded, easy to be changed)'])
end
%
%Check out that all labels are present, either thrown an error. create
%legend text.
legend_text={};n_bins=zeros(1,n_labels);
for i=min_labels:m_labels
    if all(i.*ones(l_labels,1)-labels)
        error('Labels are not natural number increasing by "1". Please create consecutive labels in previous-to-last column of input data'),
    end
    legend_text=[legend_text,['Task-phase ',num2str(i)]];
    n_bins(i)=length(labels==i);
end
%Pauses in the video decrease from 0.05, linearly
pauses=(0.02*min(n_bins).^20)./n_bins.^20;

%Centering and normalizing data for better visualization
[n,m]=size(data);
if ~(n==l_labels), error('Labels to be displayed does not match the number of patterns'), end
minimum=min(min(data));
maximum=max(max(data));
total_max=abs(minimum)+maximum;
centered=(data+abs(minimum).*ones(n,m))./(total_max);
i=1; count=1;x_prev=centered(1,1);y_prev=centered(1,2);z_prev=centered(1,3);
%
i=1; count=1;x_prev=centered(1,1);y_prev=centered(1,2);z_prev=centered(1,3);
%
max_consec_points=1;%Control of consecutive ploints re-plotted, such that plots do not take too long
%
%***Commands belonging to the real video file "video_trial_splitted_esb.m" in
%amru3_5b***
%fig=figure
%aviobj=avifile('Video_S1_Cinepak.avi','compression','Cinepak','quality',70,'fps',18);
%frames=0;
%***
%
%Centering and normalizing data for better visualization
[n,m]=size(data);
if ~(n==l_labels), error('Labels to be displayed does not match the number of patterns'), end
minimum=min(min(data));
maximum=max(max(data));
total_max=abs(minimum)+maximum;
centered=(data+abs(minimum).*ones(n,m))./(total_max);
i=1; count=1;x_prev=centered(1,1);y_prev=centered(1,2);z_prev=centered(1,3);
%
%Next variables are just for plotting the legend. Please see code below.
low_lim=-1.*ones(1,3); up_lim=ones(1,3);
for i=1:3,
    low_lim(i)=min(centered(:,i));
    up_lim(i)=max(centered(:,i));
end
step=(up_lim(3)-low_lim(3))/n_labels;
margin_1=0.1*(up_lim(1)-low_lim(1));
up_lim(1)=up_lim(1)-margin_1;%Leaving space for z-axis of the rigth plot
flags=ones(1,n_labels);
%
while count<n,
    %    frames=frames+1;
    
    if labels(count)==1
        x=[x_prev];y=[y_prev];z=[z_prev];
        consec=0;
        while (i<=n)&&(labels(i)==1),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if consec>max_consec_points,
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            subplot(2,2,index_plot)
            plot3(x_p,y_p,z_p,'Color','b','Line','-','LineWidth',1,'Marker','o','MarkerSize',1.5);%b
            hold on
            pause(pauses(1));
            %For a real video recording (not tested in here)
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig));
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        if flags(1),
            set(gca,'FontName','Arial','FontSize',8);
            tl=title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]);
            set(tl,'FontName','Arial','FontSize',10);
            xla=xlabel('PC 1');yla=ylabel('PC 2');zla=zlabel('PC 3');
            set(xla,'FontName','Arial','FontSize',10);
            set(yla,'FontName','Arial','FontSize',10);
            set(zla,'FontName','Arial','FontSize',10);
            text(up_lim(1),up_lim(2),up_lim(3),'Phase 1','Color','b','FontName','Arial','FontSize',8);
            disp('     Animation of task-phase 1...');
            flags(1)=0;
        end
        hold on
    elseif labels(count)==2,
        if flags(2),
            set(gca,'FontName','Arial','FontSize',8);
            tl=title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]);
            set(tl,'FontName','Arial','FontSize',10);
            xla=xlabel('PC 1');yla=ylabel('PC 2');zla=zlabel('PC 3');
            set(xla,'FontName','Arial','FontSize',10);
            set(yla,'FontName','Arial','FontSize',10);
            set(zla,'FontName','Arial','FontSize',10);
            text(up_lim(1),up_lim(2),up_lim(3)-step,'Phase 2','Color','r','FontName','Arial','FontSize',8);
            disp('     Animation of task-phase 2...');
            flags(2)=0;
        end
        x=[x_prev];y=[y_prev];z=[z_prev];consec=0;
        while (i<=n)&&(labels(i)==2),
            
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if consec>max_consec_points,
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            subplot(2,2,index_plot)
            plot3(x_p,y_p,z_p,'Color','r','Line','-','LineWidth',1,'Marker','o','MarkerSize',1.5);%k
            pause(pauses(2));
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig));
        end
        % reward_plotted=reward_plotted+1;
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
    elseif (labels(count)==3)&&(n_labels>2),
        if flags(3),
            set(gca,'FontName','Arial','FontSize',8);
            tl=title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]);
            set(tl,'FontName','Arial','FontSize',10);
            xla=xlabel('PC 1');yla=ylabel('PC 2');zla=zlabel('PC 3');
            set(xla,'FontName','Arial','FontSize',10);
            set(yla,'FontName','Arial','FontSize',10);
            set(zla,'FontName','Arial','FontSize',10);
            text(up_lim(1),up_lim(2),up_lim(3)-2*step,'Phase 3','Color','g','FontName','Arial','FontSize',8);
            disp('     Animation of task-phase 3...');
            flags(3)=0;
        end
        x=[x_prev];y=[y_prev];z=[z_prev];consec=0;
        while (i<=n)&&(labels(i)==3),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if (consec>max_consec_points),
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            subplot(2,2,index_plot)
            plot3(x_p,y_p,z_p,'Color','g','Line','-','LineWidth',1,'Marker','o','MarkerSize',1.5);
            pause(pauses(3));
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig))
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
    elseif (labels(count)==4)&&(n_labels>3),
        if flags(4),
            set(gca,'FontName','Arial','FontSize',8);
            tl=title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]);
            set(tl,'FontName','Arial','FontSize',10);
            xla=xlabel('PC 1');yla=ylabel('PC 2');zla=zlabel('PC 3');
            set(xla,'FontName','Arial','FontSize',10);
            set(yla,'FontName','Arial','FontSize',10);
            set(zla,'FontName','Arial','FontSize',10);
            text(up_lim(1),up_lim(2),up_lim(3)-3*step,'Phase 4','Color','k','FontName','Arial','FontSize',8);
            disp('     Animation of task-phase 4...');
            flags(4)=0;
        end
        x=[x_prev];y=[y_prev];z=[z_prev];consec=0;
        while (i<=n)&&(labels(i)==4),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if (consec>max_consec_points),
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            subplot(2,2,index_plot)
            plot3(x_p,y_p,z_p,'Color','k','Line','-','LineWidth',1,'Marker','o','MarkerSize',1.5);
            pause(pauses(4));
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig));
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
    elseif (labels(count)==5)&&(n_labels>4),
        if flags(5),
            set(gca,'FontName','Arial','FontSize',8);
            tl=title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]);
            set(tl,'FontName','Arial','FontSize',10);
            xla=xlabel('PC 1');yla=ylabel('PC 2');zla=zlabel('PC 3');
            set(xla,'FontName','Arial','FontSize',10);
            set(yla,'FontName','Arial','FontSize',10);
            set(zla,'FontName','Arial','FontSize',10);
            text(up_lim(1),up_lim(2),up_lim(3)-4*step,'Phase 4','Color','y','FontName','Arial','FontSize',8);
            disp('     Animation of task-phase 5...');
            flags(5)=0;
        end
        x=[x_prev];y=[y_prev];z=[z_prev];consec=0;
        while (i<=n)&&(labels(i)==5),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if (consec>max_consec_points),
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            subplot(2,2,index_plot)
            plot3(x_p,y_p,z_p,'Color','c','Line','-','LineWidth',1,'Marker','o','MarkerSize',1.5);
            pause(pauses(5));
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig));
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        
    elseif (labels(count)==6)&&(n_labels>5),
        if flags(6),
            set(gca,'FontName','Arial','FontSize',8);
            tl=title(['Expanded space of order ',num2str(order),' for trial #',num2str(trial)]);
            set(tl,'FontName','Arial','FontSize',10);
            xla=xlabel('PC 1');yla=ylabel('PC 2');zla=zlabel('PC 3');
            set(xla,'FontName','Arial','FontSize',10);
            set(yla,'FontName','Arial','FontSize',10);
            set(zla,'FontName','Arial','FontSize',10);
            text(up_lim(1),up_lim(2),up_lim(3)-5*step,'Phase 6','Color','c','FontName','Arial','FontSize',8);
            disp('     Animation of task-phase 6...');
            flags(6)=0;
        end
        x=[x_prev];y=[y_prev];z=[z_prev];consec=0;
        while (i<=n)&&(labels(i)==6),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if (consec>max_consec_points),
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            subplot(2,2,index_plot)
            plot3(x_p,y_p,z_p,'Color','m','Line','-','LineWidth',1,'Marker','o','MarkerSize',1.5);
            pause(pauses(1));
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig));
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
        
    elseif labels(count)==-1
        x=[x_prev];y=[y_prev];z=[z_prev];consec=0;
        while (i<=n)&&(labels(i)==-1),
            x=[x;centered(i,1)];y=[y;centered(i,2)];z=[z;centered(i,3)];
            i=i+1;consec=consec+1;
            if consec>max_consec_points,
                x_p=x(end-max_consec_points:end);
                y_p=y(end-max_consec_points:end);
                z_p=z(end-max_consec_points:end);
            else
                x_p=x;
                y_p=y;
                z_p=z;
            end
            
            plot3(x_p,y_p,z_p,'Color',[0.8,0.8,0.8],'Line','--','LineWidth',0.5,'Marker','o','MarkerSize',0.5);
            pause(pauses(6));
            %Frame(frames)=getframe;
            %aviobj = addframe(aviobj,getframe(fig));
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        hold on
    else
        error(['Task-phase number not recognized. It has to be either -1 (no any phase) or ',...
            'a natural number from 1 to ',num2str(m_labels)]),
    end
end

%Turning it around for a bit
steps=110;
angle0=get(gca,'View');
az0=angle0(1);el0=angle0(2);
for i=1:steps
    az=110*(i/steps)+az0;
    set(gca,'View',[az,el0]);
    pause(0.001);
    %if i==steps-5; pause(5); end
    %Frame(frames)=getframe;
    %aviobj = addframe(aviobj,getframe(fig));
end
steps=18;
for i=1:steps
    el=18*(i/steps)-el0;
    set(gca,'View',[az,el]);
    pause(0.001);
    %if i==steps-5; pause(5); end
    %Frame(frames)=getframe;
    %aviobj = addframe(aviobj,getframe(fig));
end
%leg=legend(legend_text);
%set(leg,'FontName','Arial','FontSize',8);
%Frame=Frame(1:length(Frame));
%aviobj=close(aviobj);
%movie(Frame)
