function display_flow2(data,labels,trial_range,order,index_plot)
%mar 011
%WARNING: Auxiliary, private file, which displays up to four
%groups flow in three dimensions (e.g. the most discriminating ones or
%orthogonalised versions of them)
%
%For clarity, only a maximum of eigth groups will be displayed in three axes.
max_groups_trajec_display=8;
scale=1;%Srinkage of the derivatived module (arrows) plotted using "quver3.m"
%        for scale=1, arrows are automatically scaled such that they do not
%       overlap and then unmmodified. scale=0 does not tmodify original arrows
%       (not advisale) and other values compress or expand all vectors afther such
%       automatical scaling.
%
m_labels=max(labels);l_labels=length(labels);min_labels=min(labels);
n_labels=m_labels-min_labels+1;
%
n=length(data(:,1));%Number of patterns
%
if n_labels>max_groups_trajec_display,
    warning('display_flow:maxEvents',['Maximum number of events displayed will be ',num2str(max_groups_trajec_display),' (hard-coded, easy to be changed)'])
end
%
%Check out that all labels are present, either thrown an error. create
%legend text.
legend_text={};
for i=min_labels:m_labels
    if all(i.*ones(l_labels,1)-labels)
        error('Labels are not natural number increasing by "1". Please create consecutive labels in previous-to-last column of input data'),
    end
end
if ~(n==l_labels), error('Labels to be displayed does not match the number of patterns'), end
%
%Normalization before plotting: +1 will be the maximum for each axis
%Small loop, no need to use 'repmat' for speed.
%Note that just the three first dimensions of the maximum discriminating subspace
%will be used, the rest are irrelevant
normalized=zeros(n,3);derivative=zeros(n-1,3);
for i=1:3,
    d=data(:,i); d_d=diff(data(:,i));
    %m_d=mean(d);m_dd=mean(d_d);
    m_d=0;m_dd=0;
    normalized(:,i)=(d-m_d)./std(d);%max(abs(d-m_d));%zscore(d);
    derivative(:,i)=(d_d-m_dd)./std(d_d);%max(abs(d_d-m_dd));%zscore(d_d);%
end
labels=labels(2:end);normalized=normalized(2:end,:);n=length(labels);
i=1; count=1;
x_prev=normalized(1,1);y_prev=normalized(1,2);z_prev=normalized(1,3);
x_d_prev=derivative(1,1);y_d_prev=derivative(1,2);z_d_prev=derivative(1,3);
%
%Next variables are just for plotting the legend. Please see code below.
low_lim=-1.*ones(1,3); up_lim=ones(1,3);
margin_1=0.2*(up_lim(1)-low_lim(1));
up_lim(1)=up_lim(1)-margin_1;%Leaving space for z-axis of the rigth plot
for i=1:3,
    [val,ind]=min(normalized(:,i));
    low_lim(i)=val+derivative(ind,i);
    [val,ind]=max(normalized(:,i));
    up_lim(i)=val+derivative(ind,i);
end
step=(up_lim(3)-low_lim(3))/n_labels;
flags=ones(1,n_labels);
while count<n,
    if labels(count)==1
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==1),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color',[0,0,.7]);
        plot3(x,y,z,'Color',[.8,.8,1],'Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(1),
            title(['Expanded space of order ',num2str(order),' for trials ',...
                num2str(trial_range(1)),' to ',num2str(trial_range(2))]),
            xlabel('DC 1');ylabel('DC 2');zlabel('DC 3');
            text(up_lim(1),low_lim(2)/2,up_lim(3),'Epoch 1','Color',[0,0,.7],'FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 1...'),
            flags(1)=0;
        end
    elseif labels(count)==2,
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==2),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color',[.1,.1,.1]);
        plot3(x,y,z,'Color',[.8,.8,.8],'Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(2),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-step,'Epoch 2','Color',[.2,.2,.2],'FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 2...'),
            flags(2)=0;
        end
    elseif (labels(count)==3)&&(n_labels>2),
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==3),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color',[.6,.6,.6]);
        plot3(x,y,z,'Color',[.9,.9,.9],'Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(3),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-2*step,'Epoch 3','Color',[.6,.6,.6],'FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 3...'),
            flags(3)=0;
        end
    elseif (labels(count)==4)&&(n_labels>3),
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==4),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color',[0,.7,0]);
        plot3(x,y,z,'Color',[0.8,1,.8],'Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(4),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-3*step,'Epoch 4','Color',[.7,0,0],'FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 4...'),
            flags(4)=0;
        end
    elseif (labels(count)==5)&&(n_labels>4),
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==5),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color',[0.7,0,0]);
        plot3(x,y,z,'Color',[1,.8,.8],'Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(5),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-4*step,'Epoch 5','Color',[.7,0,0],'FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 5...'),
            flags(5)=0;
        end
    elseif (labels(count)==6)&&(n_labels>5),
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==6),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color',[.7,.7,0]);
        plot3(x,y,z,'Color',[1,1,.8],'Line','-','LineWidth',.1);
        
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(6),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-5*step,'Epoch 6','Color',[.7,.7,0],'FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 6...'),
            flags(6)=0;
        end
    elseif (labels(count)==7)&&(n_labels>6),
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==7),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color','m');
        plot3(x,y,z,'Color','m','Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(7),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-6*step,'Epoch 7','Color','m','FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 7...'),
            flags(7)=0;
        end
    elseif (labels(count)==8)&&(n_labels>7),
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==8),
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
            i=i+1;
        end
        count=i;
        subplot(2,2,index_plot)
        quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
            'Color','c');
        plot3(x,y,z,'Color','c','Line','-','LineWidth',.1);
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
        hold on
        if flags(8),
            text(up_lim(1),low_lim(2)/2,up_lim(3)-7*step,'Epoch 8','Color','c','FontName','Arial','FontSize',8);
            disp('     Flow of task-epoch 8...'),
            flags(8)=0;
        end
    elseif labels(count)==-1
        x=x_prev;y=y_prev;z=z_prev;
        x_d=x_d_prev;y_d=y_d_prev;z_d=z_d_prev;
        while (i<=n)&&(labels(i)==-1),
            i=i+1;
            x=[x;normalized(i,1)];y=[y;normalized(i,2)];z=[z;normalized(i,3)];
            x_d=[x_d;derivative(i,1)];y_d=[y_d;derivative(i,2)];z_d=[z_d;derivative(i,3)];
        end
        %
        %Plot in ligth gray, but not during delay period
        if all(labels(count:i,end-1)~=2),
            subplot(2,2,index_plot)
            quiver3(x,y,z,x_d,y_d,z_d,scale,'Marker','o','MarkerSize',1.5,...
                'Color',[0.8,0.8,0.8]);
            plot3(x,y,z,'Color',[.8,.8,.8],'Line','-','LineWidth',.1);
            hold on
        end
        count=i;
        x_prev=x(end);y_prev=y(end);z_prev=z(end);
        x_d_prev=x_d(end);y_d_prev=y_d(end);z_d_prev=z_d(end);
    else
        error(['Task-phase number not recognized. It has to be either -1 (no any phase) or ',...
            'a natural number from 1 to ',num2str(m_labels)]),
    end
end
set(gca,'FontName','Arial','FontSize',8);
xlim([low_lim(1) up_lim(1)]);ylim([low_lim(2) up_lim(2)]);zlim([low_lim(3) up_lim(3)]);
hold off
