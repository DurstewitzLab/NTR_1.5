function display_orders(global_multiv_CE,global_multiv_DIV,twoClass_CE,twoClass_DIV,array_vals,name)
%
%This code is not commercial and thus is not guaranteed but is used on regular basis in 
%multiple datasets.
%Based on: Balaguer-Ballester E, Lapish C, Seamans JK and Durstewitz D 2011
% "Attracting Dynamics of Frontal Cortex Populations during Memory Guided Decision-Making".
% PLoS Comput. Biol.
%Distributed accourding to the General Public LIcense GNU 3.0
%http://www.gnu.org/copyleft/gpl.html
%
%by Emili Balaguer-Ballester.
%Contact adresses: eb-ballester@bournemouth.ac.uk
%                  daniel.durstewitz@zi-mannheim.de
%                  emili.balaguer@uv.es
%
%
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%
if nargin<5,
    warning('display_orders:no_orders','Eigth orders displayed'),
    array_vals=[1:8]; name='My config';
end
%    
%III.1 Multivariate
figure;
subplot(221)
errorbar(array_vals',global_multiv_CE(:,1),global_multiv_CE(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
hold on
errorbar(array_vals',global_multiv_DIV(:,1),global_multiv_DIV(:,2),'Color',[0.6,0.1,0.1],'Marker','o');
xl=xlabel('Expansion order (O)');
yl=ylabel('Segregation error (%)');
tit=title('Multivariate FDA');
leg=legend('Classification error (CE) (%)','Divergent trajectories (DIV) (%)');
set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
%
set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
%Setting the minimum/maximum for a better visualization
maCE=max(global_multiv_CE(:,1)+global_multiv_CE(:,2));maDIV=max(global_multiv_DIV(:,1)+global_multiv_DIV(:,2));
top=max(maCE,maDIV);
miCE=min(global_multiv_CE(:,1)-global_multiv_CE(:,2));miDIV=min(global_multiv_DIV(:,1)-global_multiv_DIV(:,2));
bottom=min(miCE,miDIV);
set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
hold off
%
%III.2 Two-class
subplot(222)
twoCE=mean(twoClass_CE,2);twoDIV=mean(twoClass_DIV,2);
twoCEeb=std(twoClass_CE,0,2)/sqrt(length(twoClass_CE(1,:)));
twoDIVeb=std(twoClass_DIV,0,2)/sqrt(length(twoClass_DIV(1,:)));
%
errorbar(array_vals',twoCE,twoCEeb,'Color',[0.1,0.1,0.6],'Marker','s');
hold on
errorbar(array_vals',twoDIV,twoDIVeb,'Color',[0.6,0.1,0.1],'Marker','s');
xl=xlabel('Expansion order (O)');
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
set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
hold off
%
%__________________________________
%III.3 Same but fixed span for comparison
bottom=8; top=25;%Fixed
%
%Multiv
subplot(223)
errorbar(array_vals',global_multiv_CE(:,1),global_multiv_CE(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
hold on
errorbar(array_vals',global_multiv_DIV(:,1),global_multiv_DIV(:,2),'Color',[0.6,0.1,0.1],'Marker','o');
xl=xlabel('Expansion order (O)');
yl=ylabel('Segregation error (%)');
%leg=legend('Classification error (CE) (%)','Divergent trajectories (DIV) (%)');
%set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
set(gcf,'name',name,'color',[1 1 1]);
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
hold off
%
%Two-class
subplot(224)
twoCE=mean(twoClass_CE,2);twoDIV=mean(twoClass_DIV,2);
twoCEeb=std(twoClass_CE,0,2)/sqrt(length(twoClass_CE(1,:)));
twoDIVeb=std(twoClass_DIV,0,2)/sqrt(length(twoClass_DIV(1,:)));
errorbar(array_vals',twoCE,twoCEeb,'Color',[0.1,0.1,0.5],'Marker','s');
hold on
errorbar(array_vals',twoDIV,twoDIVeb,'Color',[0.5,0.1,0.1],'Marker','s');
xl=xlabel('Expansion order (O)');
yl=ylabel('Segregation error (%)');
%leg=legend('Classification error (CE) (%)','Divergent trajectories (DIV) (%)');
%set(leg,'FontName','Helvetica','FontSize',8,'Box','off','Color','none');
set(xl,'FontName','Helvetica','FontSize',10);
set(yl,'FontName','Helvetica','FontSize',10);
set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
hold off
