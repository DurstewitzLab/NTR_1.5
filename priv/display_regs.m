function display_regs(global_multiv_CE,global_multiv_DIV,global_multiv_JS,global_multiv_L,global_multiv_L2,global_multiv_MH,...
     global_2_CE,global_2_DIV,global_2_JS,global_2_L,global_2_MH,  global_2_L2,array_vals,name)
%
if nargin<5,
    warning('display_orders:no_orders','Eigth orders displayed'),
    array_vals=(1:8); name='My config';
end
%    
%III.1 Multivariate
figure;
%
subplot(221)
errorbar(array_vals',global_multiv_CE(:,1),global_multiv_CE(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
hold on
errorbar(array_vals',global_multiv_DIV(:,1),global_multiv_DIV(:,2),'Color',[0.6,0.1,0.1],'Marker','s','LineStyle',':');
%Setting the minimum/maximum for a better visualization
maCE=max(global_multiv_CE(:,1)+global_multiv_CE(:,2));maDIV=max(global_multiv_DIV(:,1)+global_multiv_DIV(:,2));
top=max(maCE,maDIV);
miCE=min(global_multiv_CE(:,1)-global_multiv_CE(:,2));miDIV=min(global_multiv_DIV(:,1)-global_multiv_DIV(:,2));
bottom=min(miCE,miDIV);
%set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
hold off
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Segregation error (%)');
tit=title('Multivariate FDA');
leg=legend('Classification error (CE) (%)','Divergent trajectories (DIV) (%)');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%
subplot(222)
errorbar(array_vals',global_multiv_JS(:,1),global_multiv_JS(:,2),'Color',[0.6,0.1,0.1],'Marker','o');
leg=legend('Jensen-Shannon');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('JS');
tit=title('Multivariate FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%
subplot(234)
errorbar(array_vals',global_multiv_L(:,1),global_multiv_L(:,2),'Color',[0.1,0.6,0.1],'Marker','s','LineStyle','-');
leg=legend('Summed wilks lambda across significant dimensions');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Summed lambda statistic');
tit=title('Multivariate FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%
subplot(235)
errorbar(array_vals',global_multiv_L2(:,1),global_multiv_L2(:,2),'Color',[0.1,0.6,0.1],'Marker','o','LineStyle',':');
leg=legend(' Minimum wilks lambda / number of significant dim.');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Min. lambda statistic');
tit=title('Multivariate FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%
subplot(236)
errorbar(array_vals',global_multiv_MH(:,1),global_multiv_MH(:,2),'Color',[0.6,0.1,0.1],'Marker','o');
leg=legend('Averaged across-means z-distances');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Avergaged distances between groups');
tit=title('Multivariate FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
%
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%
%
%________________________________________________
figure

%III.2 Two-class
subplot(221)
%
%ce, div-2
%This is done like this just for convenience
twoCE=global_2_CE(:,1);
twoDIV=global_2_DIV(:,1);
twoCEeb=global_2_CE(:,2);
twoDIVeb=global_2_DIV(:,2);
%
errorbar(array_vals',twoCE,twoCEeb,'Color',[0.1,0.1,0.6],'Marker','o');
hold on
errorbar(array_vals',twoDIV,twoDIVeb,'Color',[0.6,0.1,0.1],'Marker','s','LineStyle',':');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Segregation error (%)');
tit=title('Pairwise FDA');
leg=legend('Classification error (CE) (%)','Divergent trajectories (DIV) (%)');
set(leg,'FontName','Helvetica','FontSize',8,'Box','off','Color','none');
%
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%Setting the minimum/maximum for a better visualization
maCE=max(twoCE+twoCEeb);maDIV=max(twoDIV+twoDIVeb);
top=max(maCE,maDIV);
miCE=min(twoCE-twoCEeb);miDIV=min(twoDIV-twoDIVeb);
bottom=min(miCE,miDIV);
%set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
hold off

%js-2
subplot(222)
errorbar(array_vals',global_2_JS(:,1),global_2_JS(:,2),'Color',[0.6,0.1,0.1],'Marker','o');
leg=legend('Pairwise Jensen-Shannon');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('JS');
tit=title('Pairwise FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);

hold on
%L-2
subplot(234)
errorbar(array_vals',global_2_L(:,1),global_2_L(:,2),'Color',[0.1,0.6,0.1],'Marker','s','LineStyle','-');
leg=legend('Summed pairwise wilks lambda across significant dimensions');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Summed lambda statistic');
tit=title('Pairwise  FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);

%L2-2
subplot(235)
errorbar(array_vals',global_2_L2(:,1),global_2_L2(:,2),'Color',[0.1,0.6,0.1],'Marker','o','LineStyle',':');
leg=legend(' Minimum pairwise wilks lambda / number of significant dim.');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Min. lambda statistic');
tit=title('Pairwise  FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%
%MH-2
subplot(236)
errorbar(array_vals',global_2_MH(:,1),global_2_MH(:,2),'Color',[0.6,0.1,0.1],'Marker','o');
leg=legend('Averaged pairwise across-means z-distances');
xl=xlabel('Reg (% mean Kernel matrix)');
yl=ylabel('Avergaged distances between groups');
tit=title('Pairwise FDA');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
%
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);

