function save_participant_metric_figure(partTbl, outDir, fname)
% save_participant_metric_figure  Bar plot of per-participant wRMSE.

    if isempty(partTbl) || height(partTbl) == 0
        return;
    end

    f = figure('Visible','off','Color','w','Name','Participant wRMSE');
    bar(partTbl.wRMSE);
    grid on;
    xlabel('participant index (participants(i))');
    ylabel('wRMSE');
    title('Per-participant weighted RMSE (unsorted)');
    exportgraphics(f, fullfile(outDir, fname), 'Resolution', 200);
    close(f);
end
