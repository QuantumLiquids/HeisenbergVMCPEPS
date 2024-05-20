
PEPSCorr;
%PEPSCorr2;
DMRGCorr;
ylim([1e-4,1]);

xticks([2,4,8,12]);
l = legend([ positive_h, h1], ...
    'PEPS, $D=8$', ...
    'DMRG, $D=15000$');
set(l,'Box','off');
set(l,'Interpreter','latex');
set(l,'Fontsize',18);
set(l,'Location','Best');