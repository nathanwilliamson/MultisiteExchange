function [] = plot_fit_and_residuals(b, I, fit)

COLORS = get(groot, 'DefaultAxesColorOrder');

fig = figure();
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.Position = [0 0 12 9];
fig.PaperPosition = [0 0 12 9];

ax1 = subplot('Position',[0.15 0.45 0.8 0.50]);
ax1.FontSize = 12;
ax1.Box = 'on';
ax1.YLabel.String = 'Intensity';
ax1.YScale = 'log';
ax1.XTickLabel = {};
ax1.YLabel.Units = 'centimeters';
ax1.YLim = [max(1e-4, 0.5*min(I)) 1.25*max(I)];

hl = line(b, I);
hl.Marker = 'o';
hl.Color = [0 0 0];
hl.LineStyle = 'none';

hl = line(b, fit.Imodel);
hl.Color = COLORS(1, :);

ax2 = subplot('Position',[0.15 0.20 0.8 0.15]);
ax2.FontSize = 12;
ax2.Box = 'on';
ax2.XLabel.String = 'b (sm^{-2})';
ax2.YLabel.String = 'Residual';
ax2.YLabel.Units = 'centimeters';
ax2.YLabel.Position(1) = ax1.YLabel.Position(1);

hl = line(b,zeros(size(b)));
hl.Color = [0 0 0];

hl = line(b, fit.residuals);
hl.Color = COLORS(1,:);
hl.Marker = 'o';
hl.LineStyle = 'none';

end