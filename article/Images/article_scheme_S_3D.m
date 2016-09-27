set(0,'defaultaxesfontsize',16,'defaulttextfontsize',16) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',5) %latex axis labels

% Bratley d=6 S_2 first order
sigma2 = 164143/2985984;
indice = 0.1791303924*sigma2;
S = @(I3, I4) min(indice./max(I3-I4.^2, eps), 1);

mu = -21/64;
I = [(sigma2 + (mu)^2) -mu];
err = [.05 .1]; 
Ireal = I - .35*err.*[1 -1];
Smax = S(I(1) - err(1), I(2) + err(2));
Smin = S(I(1) + err(1), I(2) - err(2));


step = 0.0025;
width = 0.02;
[X,Y] = meshgrid(I(1) - err(1) - width:step:I(1) + err(1) + width,...
    I(2) - err(2) - width:step:I(2) + err(2) + width);
mesh(X, Y, S(X, Y))
axis tight
xmin = min(min(X)); xmax = max(max(X));
ymin = min(min(Y)); ymax = max(max(Y));
zmin = min(min(S(X, Y))); zmax = max(max(S(X, Y)));

ax = gca;
set(ax,'TickLabelInterpreter', 'Latex');
ax.XTick = sort([I(1)-err(1) Ireal(1) I(1) I(1)+err(1)]);
ax.YTick = sort([I(2)-err(2) Ireal(2) I(2) I(2)+err(2)]);
if S(Ireal(1),Ireal(2)) < Smax & S(Ireal(1),Ireal(2)) > Smin % Not the...
    % same value as Smax and Smin for the Zaxis
    ax.ZTick = sort([Smin S(Ireal(1),Ireal(2)) (Smin+Smax)/2 S(I(1),I(2)) Smax]);
else
    ax.ZTick = sort([Smin (Smin+Smax)/2 S(I(1),I(2)) Smax]);
end
if Ireal(1) < I(1) % Ireal is indeed I, and variable I corresponds to I hat
    ax.XTickLabel = {'$\hat{I_3}-\varepsilon_{I_3}$','$I_3$','$\hat{I_3}$','$\hat{I_3}+\varepsilon_{I_3}$'};
else
    ax.XTickLabel = {'$\hat{I_3}-\varepsilon_{I_3}$','$\hat{I_3}$','$I_3$','$\hat{I_3}+\varepsilon_{I_3}$'};
end
if Ireal(2) < I(2)
    ax.YTickLabel = {'$\hat{I_4}-\varepsilon_{I_4}$','$I_4$','$\hat{I_4}$','$\hat{I_4}+\varepsilon_{I_4}$'};
else
    ax.YTickLabel = {'$\hat{I_4}-\varepsilon_{I_4}$','$\hat{I_4}$','$I_4$','$\hat{I_4}+\varepsilon_{I_4}$'};
end
if S(Ireal(1),Ireal(2)) < Smax & S(Ireal(1),Ireal(2)) > Smin % Not the...
    % same value as Smax and Smin for the Zaxis
    if S(Ireal(1),Ireal(2)) < S(I(1),I(2))
        ax.ZTickLabel = {'$S_{\min}$', '$S(I)$','$(S_{\min}+S_{\max})/2$', '$S(\hat{I})$', '$S_{\max}$'};
    else
        ax.ZTickLabel = {'$S_{\min}$', '$S(\hat{I})$', '$(S_{\min}+S_{\max})/2$', '$S(I)$', '$S_{\max}$'};
    end
else
    if S(Ireal(1),Ireal(2)) < S(I(1),I(2))
        ax.ZTickLabel = {'$S_{\min}=S(I)$', '$(S_{\min}+S_{\max})/2$', '$S(\hat{I})$', '$S_{\max}$'};
    else
        ax.ZTickLabel = {'$S_{\min}$', '$S(\hat{I})$', '$(S_{\min}+S_{\max})/2$', '$S_{\max}=S(I)$'};
    end
end
hold on
%% Vertical and horizontal lines
n = 1000; % Number of points for plotting lines
% xx = ones(n,1)*(I(1)-err(1)); % Vertical Smin corner
% yy = ones(n,1)*(I(2)-err(2));
% zz = linspace(zmin, Smax, n)'; %S(I(1)-err(1), I(2)-err(2))
% plot3(xx, yy, zz, 'k.', 'Linewidth', .5)

% xx = ones(n,1)*(I(1)+err(1)); % Vertical Smin
% yy = ones(n,1)*(I(2)-err(2));
% zz = linspace(zmin, Smin, n)'; %S(I(1)-err(1), I(2)-err(2))
% plot3(xx, yy, zz, 'k.', 'Linewidth', .5)
% 
% xx = ones(n,1)*(I(1)-err(1)); % Vertical Smax
% yy = ones(n,1)*(I(2)+err(2));
% zz = linspace(zmin, Smax, n)'; %S(I(1)-err(1), I(2)-err(2))
% plot3(xx, yy, zz, 'k.', 'Linewidth', .5)

% xx = ones(n,1)*(I(1)+err(1)); % Vertical Smax corner
% yy = ones(n,1)*(I(2)+err(2));
% zz = linspace(zmin, Smax, n)'; %S(I(1)-err(1), I(2)-err(2))
% plot3(xx, yy, zz, 'k.', 'Linewidth', .5)

% xx = linspace(I(1)-err(1), xmax, n)'; % Horizontal Smax
% yy = ones(n,1)*(I(2)+err(2));
% zz = ones(n,1)*(Smax); %S(I(1)-err(1), I(2)-err(2))
% plot3(xx, yy, zz, 'k.', 'Linewidth', .5)
% 
% xx = linspace(I(1)+err(1), xmax, n)'; % Horizontal Smin
% yy = ones(n,1)*(I(2)-err(2));
% zz = ones(n,1)*(Smin); %S(I(1)-err(1), I(2)-err(2))
% plot3(xx, yy, zz, 'k.', 'Linewidth', .5)

yy = linspace(I(2)-err(2)-width, I(2)+err(2)+width, n)'; % Plotting the (Smin+Smax)/2 line
xx = indice/((Smin+Smax)/2)+yy.^2;
zz = ones(n,1)*(Smin+Smax)/2; %S(I(1)-err(1), I(2)-err(2))
plot3(xx, yy, zz, 'k.', 'Linewidth', .5)
a = I(1) + err(1)/2; b = I(2) + err(2)/2;
% text(a, b, (Smin+Smax)/2, '$(S_{\min}+S_{\max})/2$')

xx = linspace(I(1)-err(1)-width, I(1)+err(1)+width, n)'; % Axis line for (Smin+Smax)/2 line
yy = ones(n,1)*(ymax); 
zz = ones(n,1)*(Smin+Smax)/2; %S(I(1)-err(1), I(2)-err(2))
plot3(xx, yy, zz, 'k.', 'Linewidth', .5)

aux = sqrt(abs(xmin - indice/((Smin+Smax)/2)));
yy = linspace(aux, I(2)+err(2)+width, n)'; % Axis line for (Smin+Smax)/2 line
xx = ones(n,1)*(xmin);
zz = ones(n,1)*(Smin+Smax)/2; %S(I(1)-err(1), I(2)-err(2))
plot3(xx, yy, zz, 'k.', 'Linewidth', .5)

axis([xmin xmax ymin ymax zmin zmax*1.01 ])

%area(I - err:0.02:I + err, S(I - err:0.02:I + err))
%% Limiting points on the surface and surface of possible S(I) values
marksize = 15;
plot3(Ireal(1), Ireal(2), S(Ireal(1), Ireal(2)), 'k*', 'MarkerSize', marksize)
text(Ireal(1) + 0.01, Ireal(2) + 0.01, S(Ireal(1), Ireal(2)), '$S(I)$')
plot3(I(1), I(2), S(I(1), I(2)), 'k*', 'MarkerSize', marksize)
text(I(1) + 0.01, I(2) + 0.01, S(I(1), I(2)), '$S(\hat{I})$')
plot3(I(1)-err(1), I(2)+err(2), S(I(1)-err(1), I(2)+err(2)), 'k*', 'MarkerSize', marksize)
text(I(1)-err(1) + 0.005, I(2)+err(2), S(I(1)-err(1), I(2)+err(2)) + 0.03, '$S_{\max}$')
plot3(I(1)+err(1), I(2)-err(2), S(I(1)+err(1), I(2)-err(2)), 'k*', 'MarkerSize', marksize)%, 'Color', .75*[1 1 1])
text(I(1)+err(1) - 0.005, I(2)-err(2), S(I(1)+err(1), I(2)-err(2)) + 0.1, '$S_{\min}$')%, 'Color', .75*[1 1 1])

xx = linspace(I(1)-err(1), I(1)+err(1), n)'; % Corner 1
yy = ones(n,1)*(I(2)-err(2));
zz = S(xx, yy);
plot3(xx, yy, zz, 'b.', 'Linewidth', .5)

xx = linspace(I(1)-err(1), I(1)+err(1), n)'; % Corner 2
yy = ones(n,1)*(I(2)+err(2));
zz = S(xx, yy);
plot3(xx, yy, zz, 'b.', 'Linewidth', .5)

xx = ones(n,1)*(I(1)-err(1)); % Corner 3
yy = linspace(I(2)-err(2), I(2)+err(2), n)';
zz = S(xx, yy);
plot3(xx, yy, zz, 'b.', 'Linewidth', .5)

xx = ones(n,1)*(I(1)+err(1)); % Corner 3
yy = linspace(I(2)-err(2), I(2)+err(2), n)';
zz = S(xx, yy);
plot3(xx, yy, zz, 'b.', 'Linewidth', .5)

hold off

[az, e1] = view;
% az = az + 80;
% az = 40.5; e1 = 46;
az = -15.5; e1 = 30;
view(az,e1)

colormap autumn
print '-depsc2' -opengl  'estimator_3d_small.eps'