%%
folder='E:\';
f=figure;
f.PaperPosition=[1 1 7 6];
set(gca,'fontsize',9)
hold all 
step=2;
range=400;
X=1./([250 300 350]+273.15);

plot(X,[2.90E-15 1.00E-14 2.17E-14],'-o',...
    'Markersize',4,'Markerfacecolor','r',...
    'color',[0 0 0],'Markeredgecolor','k');
box on

ax=gca;
ax.XLabel.String='1/T (K^{-1})';
ax.XLabel.FontSize=10;
ax.YLabel.String='D (cm^2/s)';
ax.YLabel.FontSize=10;
xlim([1.55e-3 1.95e-3]);
ylim([2e-15 2.5e-14]);
set(gca, 'YScale', 'log');

print(f,'-dpng','-r300',[folder 'layer intensity.png']);