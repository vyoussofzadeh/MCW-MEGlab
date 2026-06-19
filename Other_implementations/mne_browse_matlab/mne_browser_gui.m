function mne_browser_gui(filename)
% MNE_BROWSER_GUI – Raw MEG/EEG browser (MATLAB UIFigure)
%
% v1.30 | 2025-06-14
%   • Fixes “colon-operator” bugs in helper creators (addField / addDrop).
%   • X-axis now always shows the real time range (e.g. 25–30 s).
%   • Same MAG / GRAD / EEG / ECG / Other scaling and Auto-match logic as v1.29.
%   • Classic label column lines up; colour cycling is MATLAB default.
% -------------------------------------------------------------------------

assert(exist('ft_read_header','file')==2,'FieldTrip (or ft_read_header) not found.');

%% ---------------------------------------------------------------------  HEADER
hdr      = ft_read_header(filename);
fs       = hdr.Fs;
labels   = hdr.label;                   % original labels (with spaces)
labNS    = strrep(labels,' ','');       % no-space copies (for matching)
nsamples = hdr.nSamples * hdr.nTrials;
T        = nsamples / fs;

[roiDefs, roiNames] = mne_default_montages(labels);

%% ---------------------------------------------------------------------  CONSTANTS
margin  = 10;
listW   = 160;
panelW  = 320;
sliderH = 30;

%% ---------------------------------------------------------------------  FIGURE + MAIN UI
fig = uifigure('Name','MNE Browser','Position',[60 80 1500 780]);
fig.AutoResizeChildren = 'off';

ui = struct();

% channel list
ui.lst = uilistbox(fig,'Items',labels,'MultiSelect','on',...
                   'Value',labels(1:min(20,end)));

% plotting axes
ui.ax       = uiaxes(fig,'Box','on','YTick',[]); hold(ui.ax,'on');
ui.selLine  = [];    % cursor
ui.curLabel = uilabel(fig,'Text','','FontWeight','bold');

% time slider (bottom)
ui.sld = uislider(fig,'Limits',[0 T],'MajorTicks',[],'MinorTicks',[]);

%% ---------------------------------------------------------------------  SETTINGS PANEL
setP  = uipanel(fig,'Title','Settings');
Y0    = 700;
dy    = 26;
X0    = 10;
Y     = Y0;

addField('Window (s):','lenEdit',5,0.5,60);
addDrop ('ROI / Montage:','roiDrop',roiNames,'Manual');
addCheck('Full view','flvChk');
addCheck('Classic labels','classicChk');

Y = Y - dy*2;
uilabel(setP,'Text','Filter','FontWeight','bold','Position',[X0 Y 120 22]);
addCheck('Enable filter','fltChk');
addField('Low (Hz):' ,'lowEdit' ,1   ,0     ,fs/2-1e-3);
addField('High (Hz):','highEdit',40  ,0.1   ,fs/2-1e-3);

Y = Y - dy*2;
uilabel(setP,'Text','Scale factors','FontWeight','bold','Position',[X0 Y 120 22]);
addField('MAG (…1):'  ,'magScaleEdit' ,1   ,0.001,1000);
addField('GRAD (…2/3):','gradScaleEdit',1  ,0.001,1000);
addField('EEG:'        ,'eegScaleEdit' ,1  ,0.001,1000);
addField('ECG:'        ,'ecgScaleEdit' ,1e-3,1e-9,1000);
addField('Other:'      ,'othScaleEdit' ,1  ,0.001,1000);

Y = Y - dy;
ui.autoBtn = uibutton(setP,'Text','Auto-match',...
    'ButtonPushedFcn',@autoScaleCb,...
    'Position',[X0 Y 110 22]);

Y = Y - dy;
uilabel(setP,'Text','Vert. Gain','FontWeight','bold','Position',[X0 Y 120 22]);
addField('× σ:','gainEdit',2.5,0.5,8);

Y = Y - dy*2;
uilabel(setP,'Text','Navigation','FontWeight','bold','Position',[X0 Y 120 22]);
addField('Go to (s):','timeEdit',0,0,T);

Y = Y - dy;
ui.prevBtn = uibutton(setP,'Text','⟵ Prev',...
    'Position',[X0      Y 80 22],'ButtonPushedFcn',@(~,~)shiftWin(-1));
ui.nextBtn = uibutton(setP,'Text','Next ⟶',...
    'Position',[X0+90  Y 80 22],'ButtonPushedFcn',@(~,~)shiftWin(+1));

%% ---------------------------------------------------------------------  LAYOUT + RESIZE
reflow();
fig.SizeChangedFcn = @(~,~)reflow();

    function reflow()
        W = fig.Position(3);  H = fig.Position(4);
        ui.lst.Position   = [margin margin+sliderH+margin listW H-3*margin-sliderH];
        setP.Position     = [W-panelW-margin margin+sliderH+margin panelW H-3*margin-sliderH];

        axX = margin + listW + margin;
        axW = W - axX - panelW - 2*margin;
        axH = H - 3*margin - sliderH;
        axY = margin + sliderH + margin;
        ui.ax.Position  = [axX axY axW axH];

        ui.sld.Position = [axX margin axW sliderH-10];
        ui.curLabel.Position = [axX margin+sliderH-15 240 22];
    end

%% ---------------------------------------------------------------------  TIME EDIT / SLIDER SYNC
ui.timeEdit.ValueChangedFcn = @(src,~)jumpTo(src.Value);
ui.sld.ValueChangedFcn      = @(src,~)jumpTo(src.Value);

%% ---------------------------------------------------------------------  MOUSE CLICK → CURSOR
try, fig.WindowButtonDownFcn = @figClick; catch, end
    function figClick(~,~)
        cp = ui.ax.CurrentPoint;
        if all(cp(1,1:2) >= [ui.ax.XLim(1) ui.ax.YLim(1)]) && ...
           all(cp(1,1:2) <= [ui.ax.XLim(2) ui.ax.YLim(2)])
            jumpTo(cp(1,1));
        end
    end

%% ---------------------------------------------------------------------  KEYBOARD SHORTCUTS
try, fig.KeyPressFcn = @keyCb; catch, end
    function keyCb(~,ev)
        if strcmp(ev.Key,'leftarrow'),  shiftWin(-1); end
        if strcmp(ev.Key,'rightarrow'), shiftWin(+1); end
    end

%% ---------------------------------------------------------------------  NAVIGATION HELPERS
    function shiftWin(sign)
        jumpTo(ui.timeEdit.Value + sign * ui.lenEdit.Value);
    end

    function jumpTo(newT)
        newT = max(0, min(newT, T-ui.lenEdit.Value));
        ui.timeEdit.Value = newT;
        ui.sld.Value      = newT;
        updatePlot();
    end

%% ---------------------------------------------------------------------  CALLBACK REGISTRATION
allCb = [ui.lst ui.sld ui.lenEdit ui.fltChk ui.lowEdit ui.highEdit ...
         ui.magScaleEdit ui.gradScaleEdit ui.eegScaleEdit ui.ecgScaleEdit ui.othScaleEdit ui.gainEdit ...
         ui.flvChk ui.roiDrop ui.classicChk];
set(allCb,'ValueChangedFcn',@(~,~)updatePlot());

updatePlot();   % first draw

%% ---------------------------------------------------------------------  AUTO-MATCH BUTTON
    function autoScaleCb(~,~)
        t0  = ui.timeEdit.Value;
        win = ui.lenEdit.Value;
        iBeg = max(1,floor(t0*fs)+1);
        iEnd = min(nsamples,floor((t0+win)*fs));

        % channel indices in current view
        chIdx = getCurrentIdx();
        if isempty(chIdx), return; end

        datRaw = ft_read_data(filename,'header',hdr,...
                              'begsample',iBeg,'endsample',iEnd,...
                              'chanindx',chIdx);
        chRMS  = std(datRaw,0,2,'omitnan');
        getMed = @(sel) median(chRMS(sel),'omitnan');

        isMag  = endsWith(labNS(chIdx),'1');
        isGrad = endsWith(labNS(chIdx),'2') | endsWith(labNS(chIdx),'3');
        mMAG   = getMed(isMag);
        mGRAD  = getMed(isGrad);
        mEEG   = getMed(startsWith(labNS(chIdx),'EEG'));
        mECG   = getMed(startsWith(labNS(chIdx),'ECG'));
        mOTH   = getMed(~(isMag | isGrad | startsWith(labNS(chIdx),'EEG') | startsWith(labNS(chIdx),'ECG')));

        ref  = median([mMAG mGRAD mEEG mECG mOTH],'omitnan'); if isnan(ref)||ref==0, ref = 1; end
        clip = @(v,lim) max(min(v,lim(2)),lim(1));

        if mMAG >0,  ui.magScaleEdit.Value  = clip(ref/mMAG , ui.magScaleEdit.Limits); end
        if mGRAD>0, ui.gradScaleEdit.Value = clip(ref/mGRAD, ui.gradScaleEdit.Limits); end
        if mEEG >0, ui.eegScaleEdit.Value  = clip(ref/mEEG , ui.eegScaleEdit.Limits); end
        if mECG >0, ui.ecgScaleEdit.Value  = clip(ref/mECG , ui.ecgScaleEdit.Limits); end
        if mOTH >0, ui.othScaleEdit.Value  = clip(ref/mOTH , ui.othScaleEdit.Limits); end

        updatePlot();
    end

%% ---------------------------------------------------------------------  PLOT UPDATE
    function updatePlot()
        reflow();

        % time window -------------------------------------------------
        t0  = ui.timeEdit.Value;
        win = ui.lenEdit.Value;
        iBeg = max(1, floor(t0*fs)+1);
        iEnd = min(nsamples, floor((t0+win)*fs));

        % channel selection ------------------------------------------
        chIdx = getCurrentIdx();
        if isempty(chIdx), cla(ui.ax); return; end

        % data read ---------------------------------------------------
        dat = ft_read_data(filename,'header',hdr,...
                           'begsample',iBeg,'endsample',iEnd,...
                           'chanindx',chIdx);

        % optional filter
        if ui.fltChk.Value
            lo = ui.lowEdit.Value;  hi = ui.highEdit.Value;
            if lo < hi && hi < fs/2
                [b,a] = butter(4,[lo hi]/(fs/2));
                dat   = filtfilt(b,a,dat')';
            end
        end

        % modality-specific scaling ----------------------------------
        for k = 1:size(dat,1)
            lab = labNS{chIdx(k)};
            scale = ui.othScaleEdit.Value;
            if endsWith(lab,'1')
                scale = ui.magScaleEdit.Value;
            elseif endsWith(lab,'2') || endsWith(lab,'3')
                scale = ui.gradScaleEdit.Value;
            elseif startsWith(lab,'EEG')
                scale = ui.eegScaleEdit.Value;
            elseif startsWith(lab,'ECG')
                scale = ui.ecgScaleEdit.Value;
            end
            dat(k,:) = dat(k,:) * scale;
        end

        % vertical separation & plot ---------------------------------
        sep  = ui.gainEdit.Value * (std(dat(:),'omitnan') + eps);
        yOff = ((size(dat,1)-1):-1:0) * sep;
        tvec = (iBeg:iEnd)/fs;        % absolute time axis

        cla(ui.ax);
        for k = 1:size(dat,1)
            plot(ui.ax,tvec,dat(k,:) + yOff(k)); hold(ui.ax,'on');
        end
        ui.ax.XLim = [t0 t0+win];

        if ui.classicChk.Value
            ui.ax.YTick      = fliplr(yOff);
            ui.ax.YTickLabel = flip(labels(chIdx));
            ui.ax.FontSize   = 10;
        else
            ui.ax.YTick      = [];
            ui.ax.YTickLabel = {};
        end
        ui.curLabel.Text = sprintf('t = %.2f – %.2f s', t0, t0+win);

        % cursor line ------------------------------------------------
        if ~isempty(ui.selLine) && isgraphics(ui.selLine), delete(ui.selLine); end
        if t0 <= ui.timeEdit.Value && ui.timeEdit.Value <= t0+win
            ui.selLine = line(ui.ax,[ui.timeEdit.Value ui.timeEdit.Value], ui.ax.YLim,...
                              'Color','k','LineWidth',1.5);
        else
            ui.selLine = [];
        end
    end

%% ---------------------------------------------------------------------  CURRENT IDX HELPER
    function chIdx = getCurrentIdx()
        if ui.flvChk.Value
            chIdx = 1:numel(labels);
            ui.lst.Enable = 'off';
        elseif ~strcmp(ui.roiDrop.Value,'Manual')
            tgtNS = strrep(roiDefs(ui.roiDrop.Value),' ','');
            chIdx = find(ismember(labNS, tgtNS));
            ui.lst.Enable = 'off';
        else
            chIdx = find(ismember(labels, ui.lst.Value));
            ui.lst.Enable = 'on';
        end
    end

%% ---------------------------------------------------------------------  HELPER-CREATOR FUNCTIONS
    function addField(lbl,name,def,minv,maxv)
        Y = Y - dy;
        uilabel(setP,'Text',lbl,'Position',[X0 Y 95 22]);
        ui.(name) = uieditfield(setP,'numeric',...
            'Value',def,'Limits',[minv maxv],...
            'Position',[115 Y 90 22]);
    end

    function addDrop(lbl,name,items,val)
        Y = Y - dy;
        uilabel(setP,'Text',lbl,'Position',[X0 Y 95 22]);
        ui.(name) = uidropdown(setP,...
            'Items',items,'Value',val,...
            'Position',[115 Y 180 22]);
    end

    function addCheck(lbl,name)
        Y = Y - dy;
        ui.(name) = uicheckbox(setP,'Text',lbl,...
            'Position',[X0 Y 200 22]);
    end
end
