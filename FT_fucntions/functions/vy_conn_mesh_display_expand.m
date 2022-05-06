function tfilenameSURF=vy_conn_mesh_display_expand(filenameSURF, idx)
tfilenameSURF=cellstr(conn_expandframe(filenameSURF));
try
    trefnames=regexp(fileread(conn_prepend('',filenameSURF,'.txt')),'\n*','split');
    trefnames=trefnames(cellfun('length',trefnames)>0);
    if numel(trefnames)>1&&numel(tfilenameSURF)==1 %3d-atlas
        tempvol=spm_vol(char(tfilenameSURF));
        tempdata=spm_read_vols(tempvol);
        maxdata=max(tempdata(:));
        if numel(trefnames)==maxdata
            idata=idx;%listdlg('liststring',trefnames,'selectionmode','multiple','initialvalue',1:numel(trefnames),'promptstring','Select ROI:','ListSize',[300 200]);
%             idata=listdlg('liststring',trefnames,'selectionmode','multiple','initialvalue',1:numel(trefnames),'promptstring','Select ROI:','ListSize',[300 200]);
            if isempty(idata), return; end
            tfilenameSURF=num2cell(struct('filename',char(filenameSURF),'data',arrayfun(@(n)tempdata==n,idata,'uni',0),'vol',tempvol));
        end
    elseif numel(trefnames)>1 && numel(trefnames)==numel(tfilenameSURF) %4d-atlas
        idata=listdlg('liststring',trefnames,'selectionmode','multiple','initialvalue',1:numel(trefnames),'promptstring','Select ROI:','ListSize',[300 200]);
        if isempty(idata), return; end
        tfilenameSURF=tfilenameSURF(idata);
    end
end
end