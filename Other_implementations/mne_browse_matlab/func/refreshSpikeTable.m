function refreshSpikeTable
if isfield(ui,'spikeWin') && isvalid(ui.spikeWin)
    tbl = findobj(ui.spikeWin,'Type','uitable');
    if ~isempty(tbl)
        tbl.Data = num2cell(ui.spikeT(:));
    end
end
end
