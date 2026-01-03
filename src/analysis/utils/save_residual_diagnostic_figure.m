function save_residual_diagnostic_figure(resTbl, outDir, fname)
% save_residual_diagnostic_figure  2x2 residual diagnostics figure.

    if isempty(resTbl) || height(resTbl) == 0
        return;
    end

    r    = resTbl.residual;
    yhat = resTbl.y_hat;
    t    = resTbl.t_s;
    kind = resTbl.kind;

    f = figure('Visible','off','Color','w','Name','Residual diagnostics');
    tiledlayout(2,2,'Padding','loose','TileSpacing','loose');

    % 1) Histogram
    nexttile;
    histogram(r, 30);
    grid on;
    xlabel('residual (y_{obs} - y_{hat})');
    ylabel('count');
    title('Residual histogram');

    % 2) Boxplot by kind (robust ordering)
    nexttile;
    desired = ["t40_post","t14_mid1","t14_mid2","probe"];
    kindCat = categorical(kind, desired, 'Ordinal', true);

    % Remove undefined (kinds not in desired)
    keep = ~isundefined(kindCat);
    boxplot(r(keep), kindCat(keep));
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 25;
    xlabel('measurement type');
    ylabel('residual');
    title('Residuals by measurement type');

    % 3) Residual vs predicted
    nexttile;
    plot(yhat, r, '.', 'MarkerSize', 8);
    grid on;
    xlabel('y_{hat}');
    ylabel('residual');
    title('Residual vs predicted');

    % 4) Residual vs time
    nexttile;
    plot(t, r, '.', 'MarkerSize', 8);
    grid on;
    xlabel('time [s]');
    ylabel('residual');
    title('Residual vs time');

    exportgraphics(f, fullfile(outDir, fname), 'Resolution', 200);
    close(f);
end
