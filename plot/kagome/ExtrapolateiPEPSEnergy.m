L = [6,8,10,12];
e_site = [-0.438356, -0.437953, -0.437580, -0.437383];
e_site_all = [-0.409008, -0.415044,-0.418851, -0.421563];
plot(1./L,e_site,'o');hold on;
plot(1./L, e_site_all,'o'); hold on;

e_dider = -0.436;

fit_x = 1./L;
p = fit(fit_x',e_site','poly1');
fprintf('e_ex=%.5f\n',p.p2);
x = 0:0.01:max(fit_x);
h1=loglog(x ,p.p1*x+p.p2,'-.');

fit_x = 1./L;
p = fit(fit_x(2:end)',e_site_all(2:end)','poly1');
fprintf('e_ex=%.5f\n',p.p2);
x = 0:0.01:max(fit_x);
h2=loglog(x ,p.p1*x+p.p2,'-.');

h3=yline(e_dider, 'b--');

l=legend([h1, h2,h3], {'persite energy in bulk', 'persite energy totally', ...
    'Didier energy for $D=8$'});
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',18);
set(l,'Location','NorthEast');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$1/L$','Interpreter','latex');
ylabel('persite energy','Interpreter','latex');
% ylabel('$S_z$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);