function spikeSelect(src,ev)
if isempty(ev.Indices), return; end
tSel = src.Data{ev.Indices(1)};
jumpTo(tSel);
delete(findall(ui.ax,'Tag','SpikeSel'));
line(ui.ax,[tSel tSel], ui.ax.YLim, ...
    'Color','r','LineWidth',1.5,'Tag','SpikeSel');
end
