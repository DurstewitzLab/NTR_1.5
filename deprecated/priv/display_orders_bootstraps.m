function display_orders_bootstraps(multiv_CE,multiv_DIV,twoclass_CE,twoclass_DIV,array_vals,name,fixed_span)
%
if nargin<6, name=''; end
if nargin<7, fixed_span=[]; end
if fixed_span
    %Setting the minimum/maximum for a better visualization
    maCE=max(multiv_CE(:,1)+multiv_CE(:,2));maDIV=max(multiv_DIV(:,1)+multiv_DIV(:,2));
    top=max(maCE,maDIV);
    miCE=min(multiv_CE(:,1)-multiv_CE(:,2));miDIV=min(multiv_DIV(:,1)-multiv_DIV(:,2));
    bottom=min(miCE,miDIV);
else
    %bottom=fixed_span(1); top=fixed_span(2);
    bottom=8; top=25;
    warning('display_orders_bootstraps:fixed_span','Span of error (y-axes) was fixed for all graphics. "fixed_span" argument controls this effect'),
end
if nargin<5,
    warning('display_orders:no_orders','Ten orders displayed'),
    array_vals=(1:10); name='My config';
end
n=length(array_vals);
%
figure;
%I. Multivariate display
%I.1 CE with bootstraps
    subplot(221)
    errorbar(array_vals',multiv_CE(:,1),multiv_CE(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
    hold on
    errorbar(array_vals',multiv_CE(:,3),multiv_CE(:,4),'Color',[0.6,0.1,0.1],'Marker','s');
    xl=xlabel('Expansion order (O)');
    yl=ylabel('Segregation error (%)');
    tit=title('Multivariate FDA');
    leg=legend('Original data',['Bootstraps (n=',num2str(n)']);
    set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
    set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
    set(gcf,'name',name,'color',[1 1 1]);
    set(xl,'FontName','Helvetica','FontSize',10);
    set(yl,'FontName','Helvetica','FontSize',10);
    set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
    hold off
%
%I.2 DIV with bootstraps
    subplot(222)
    errorbar(array_vals',multiv_DIV(:,1),multiv_DIV(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
    hold on
    errorbar(array_vals',multiv_DIV(:,3),multiv_DIV(:,4),'Color',[0.6,0.1,0.1],'Marker','s');
    xl=xlabel('Expansion order (O)');
    yl=ylabel('Divergent trajectories (%)');
    tit=title('Multivariate FDA');
    leg=legend('Original data',['Bootstraps (n=',num2str(n)']);
    set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
    set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
    set(gcf,'name',name,'color',[1 1 1]);
    set(xl,'FontName','Helvetica','FontSize',10);
    set(yl,'FontName','Helvetica','FontSize',10);
    set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
    hold off
%    
%III.2 Two-class
%I.1 CE with bootstrap
    subplot(223)
    errorbar(array_vals',twoclass_CE(:,1),twoclass_CE(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
    hold on
    errorbar(array_vals',twoclass_CE(:,3),twoclass_CE(:,4),'Color',[0.6,0.1,0.1],'Marker','s');
    xl=xlabel('Expansion order (O)');
    yl=ylabel('Segregation error (%)');
    tit=title('Pairwise FDA');
    leg=legend('Original data',['Bootstraps (n=',num2str(n)']);
    set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
    set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
    set(gcf,'name',name,'color',[1 1 1]);
    set(xl,'FontName','Helvetica','FontSize',10);
    set(yl,'FontName','Helvetica','FontSize',10);
    set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
    hold off
%I.2 DIV with bootstraps
    subplot(224)
    errorbar(array_vals',twoclass_DIV(:,1),twoclass_DIV(:,2),'Color',[0.1,0.1,0.6],'Marker','o');
    hold on
    errorbar(array_vals',twoclass_DIV(:,3),twoclass_DIV(:,4),'Color',[0.6,0.1,0.1],'Marker','s');
    xl=xlabel('Expansion order (O)');
    yl=ylabel('Divergent trajectories (%)');
     tit=title('Pairwise FDA');
    leg=legend('Original data',['Bootstraps (n=',num2str(n)']);
    set(leg,'FontName','Helvetica','FontSize',8, 'Box','off','Color','none');
    set(tit,'FontName','Helvetica','FontSize',10,'FontAngle','Italic');
    set(gcf,'name',name,'color',[1 1 1]);
    set(xl,'FontName','Helvetica','FontSize',10);
    set(yl,'FontName','Helvetica','FontSize',10);
    set(gca,'FontName','Helvetica','FontSize',8,'YLim',[bottom top],'YScale','log');
    hold off
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

