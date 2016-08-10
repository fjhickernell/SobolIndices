    set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
    set(0,'defaultLineLineWidth',3) %thick lines
    set(0,'defaultTextInterpreter','latex') %latex axis labels
    set(0,'defaultLegendInterpreter','latex') %latex axis labels
    set(0,'defaultLineMarkerSize',40) %latex axis labels

S = @(x) exp(-(x-1).^2) + 1;

I = 0.90;
err = 0.2;
Smax = S(1); xmax = 1;
Smin = S(I - err); xmin = I - err;


x = I - err - 0.06:0.02:I + err + 0.06;
ymin = min(S(x));
ymax = max(S(x));

plot(x,S(x))
ax = gca;
set(ax,'TickLabelInterpreter', 'Latex');
ax.XTick = [I-err I I+err];
ax.YTick = [Smin (Smin+Smax)/2 S(I) Smax];
ax.XTickLabel = {'$\hat{I}-err$','$\hat{I}$','$\hat{I}+err$'};
ax.YTickLabel = {'$S_{\min}$', '$(S_{\min}+S_{\max})/2$', '$S(\hat{I})$', '$S_{\max}$'};
hold on
plot([xmin xmin], [ymin Smin], 'k--')
plot([x(1) xmin], [Smin Smin], 'k--')
plot([I I], [ymin S(I)], 'k--')
plot([x(1) I], [S(I) S(I)], 'k--')
plot([xmax xmax], [ymin Smax], 'k--')
plot([x(1) xmax], [Smax Smax], 'k--')
hold off
axis([x(1) x(end) ymin ymax*1.01])