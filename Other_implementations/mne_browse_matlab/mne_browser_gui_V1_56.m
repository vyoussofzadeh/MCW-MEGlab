function mne_browser_gui_V1_56(filename)
% MNE_BROWSER_GUI – Raw MEG/EEG browser (MATLAB UIFigure)
%
% v1.56 | 2025-06-15
%   • Adds one-click “Detect spikes” (80-250 Hz, z-score > 5).
%   • Spike times persist while the session is open and are drawn as
%     magenta dashed lines on the trace.
% -------------------------------------------------------------------------

assert(exist('ft_read_header','file')==2,'FieldTrip (ft_read_header) not found.');

addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Other_implementations/mne_browse_matlab/func')

%% ----------------------------- HEADER ---------------------------------
hdr      = ft_read_header(filename);
fs       = hdr.Fs;
labels   = hdr.label;
labNS    = strrep(labels,' ','');          % no-space copies
nsamples = hdr.nSamples * hdr.nTrials;
T        = nsamples / fs;

[roiDefs, roiNames] = mne_default_montages(labels);

%% ---------------------------- CONSTANTS -------------------------------
margin  = 10;  listW = 160;  panelW = 320;  sliderH = 30;

%% ---------------------- FIGURE & MAIN CONTROLS ------------------------
fig = uifigure('Name','MNE Browser','Position',[60 80 1500 780]);
fig.AutoResizeChildren = 'off';


ui = struct();
ui.spikeT = [];                         % <-- stores spike times (s)
ui.autoScale = true;        % <-- put it here



ui.lst = uilistbox(fig,'Items',labels,'MultiSelect','on',...
    'Value',labels(1:min(end,20)));

ui.ax       = uiaxes(fig,'Box','on','YTick',[]); hold(ui.ax,'on');
ui.selLine  = [];
ui.curLabel = uilabel(fig,'Text','','FontWeight','bold');

ui.sld = uislider(fig,'Limits',[0 T],'MajorTicks',[],'MinorTicks',[]);


% ─── Spike menu in the menubar ─────────────────────────────────────
mSpike = uimenu(fig,'Text','Spike');
uimenu(mSpike,'Text','Detect spikes…', ...
    'MenuSelectedFcn',@(~,~)openSpikePanel( ...
    ui, filename, hdr, fs, nsamples, ...
    @getCurrentIdx, @updatePlot));

%% --------------------------- SETTINGS ---------------------------------
setP = uipanel(fig,'Title','Settings');
Y0  = 700;  dy = 26;  X0 = 10;  Y = Y0;

Y = Y-dy*1;  uilabel(setP,'Text','Montage','FontWeight','bold','Position',[X0 Y 120 22]);

Y = Y-dy*1; addDrop (setP, X0,Y,'ROI / Montage:','roiDrop',roiNames,'Manual');
addCheck('Full view','flvChk');
addCheck('Classic labels','classicChk');

Y = Y-dy*1;  uilabel(setP,'Text','Filter','FontWeight','bold','Position',[X0 Y 120 22]);
addCheck('Enable filter','fltChk');
addField('Low (Hz):','lowEdit',1  ,0   ,fs/2-1e-3);
addField('High (Hz):','highEdit',40 ,0.1,fs/2-1e-3);

Y = Y-dy*1;  uilabel(setP,'Text','Scale factors','FontWeight','bold','Position',[X0 Y 120 22]);
% addField('MAG (…1):'  ,'magScaleEdit' ,1   ,0.001,1000);
% addField('GRAD (…2/3):','gradScaleEdit',1  ,0.001,1000);
% addField('EEG:'        ,'eegScaleEdit' ,1  ,0.001,1000);
% addField('ECG:'        ,'ecgScaleEdit' ,1e-3,1e-9,1000);
% addField('Other:'      ,'othScaleEdit' ,1  ,0.001,1000);

% ---- Scale block -------------------------------------------------
addField('MAG (…1):'  ,'magScaleEdit' ,1,0.001,10000,true);
addField('GRAD (…2/3):','gradScaleEdit',1,0.001,10000,true);
addField('EEG:'        ,'eegScaleEdit' ,1,0.001,10000,true);
addField('ECG:'        ,'ecgScaleEdit' ,1e-3,1e-9,10000,true);
addField('Other:'      ,'othScaleEdit' ,1,0.001,10000,true);

Y = Y - dy;
ui.autoBtn = uibutton(setP, ...
    'Text','Auto-match', ...
    'Position',[X0 Y 110 22], ...
    'ButtonPushedFcn',@enableAuto);

    function enableAuto(~,~)
        ui.autoScale = true;
        autoScaleNow();
        updatePlot();
    end

% -----------------------------------------------------------------------

Y = Y-dy;  uilabel(setP,'Text','Vert. Gain','FontWeight','bold','Position',[X0 Y 120 22]);
addField('× σ:','gainEdit',7,0.5,8);

Y = Y-dy*1;  uilabel(setP,'Text','Navigation','FontWeight','bold','Position',[X0 Y 120 22]);
addField('Go to (s):','timeEdit',0,0,T);
addField('Window (s):','lenEdit',5,0.5,60);

Y = Y-dy;
ui.prevBtn = uibutton(setP,'Text','⟵ Prev','Position',[X0 Y 80 22],...
    'ButtonPushedFcn',@(~,~)shiftWin(-1));
ui.nextBtn = uibutton(setP,'Text','Next ⟶','Position',[X0+90 Y 80 22],...
    'ButtonPushedFcn',@(~,~)shiftWin(+1));

%% --------------------------- LAYOUT -----------------------------------
reflow();
fig.SizeChangedFcn = @(~,~)reflow();


%% --------------------- NAVIGATION BASICS ------------------------------
ui.timeEdit.ValueChangedFcn = @(s,~)jumpTo(s.Value);
ui.sld.ValueChangedFcn      = @(s,~)jumpTo(s.Value);
try, fig.WindowButtonDownFcn = @mouseClick; catch, end
    function mouseClick(~,~)
        cp=ui.ax.CurrentPoint;
        if cp(1,1)>=ui.ax.XLim(1)&&cp(1,1)<=ui.ax.XLim(2)
            jumpTo(cp(1,1)); end
    end
try, fig.KeyPressFcn = @(~,ev)keyHit(ev); catch, end
    function keyHit(ev)
        switch ev.Key
            case 'leftarrow',  shiftWin(-1);
            case 'rightarrow', shiftWin(+1);
        end
    end
    function shiftWin(dir), jumpTo(ui.timeEdit.Value+dir*ui.lenEdit.Value); end
    function jumpTo(t)
        t=max(0,min(t,T-ui.lenEdit.Value));
        ui.timeEdit.Value=t; ui.sld.Value=t; updatePlot();
    end

%% ----------------- CALLBACK REGISTRATION ------------------------------
% allCb = [ui.lst ui.sld ui.lenEdit ui.fltChk ui.lowEdit ui.highEdit ...
%     ui.magScaleEdit ui.gradScaleEdit ui.eegScaleEdit ui.ecgScaleEdit ui.othScaleEdit ui.gainEdit ...
%     ui.flvChk ui.roiDrop ui.classicChk];
% set(allCb,'ValueChangedFcn',@(~,~)updatePlot());

% everything *except* the five scale edits
allCb = [ui.lst ui.sld ui.lenEdit ui.fltChk ui.lowEdit ui.highEdit ...
    ui.gainEdit ui.flvChk ui.roiDrop ui.classicChk];
set(allCb,'ValueChangedFcn',@(~,~)updatePlot());

updatePlot();   % initial draw

    function addDrop(parent,x0,y,lbl,name,items,val)
        uilabel(parent,'Text',lbl,'Position',[x0 y 95 22]);
        ui.(name) = uidropdown(parent, ...
            'Items',items,'Value',val, ...
            'Position',[x0+105 y 180 22]);
    end

    function addCheck(lbl,name)
        Y=Y-dy; ui.(name)=uicheckbox(setP,'Text',lbl,...
            'Position',[X0 Y 200 22]);
    end
%     function addField(lbl,name,def,minv,maxv,isScale)
%         if nargin<6, isScale = false; end          % default
%         Y = Y - dy;
%         uilabel(setP,'Text',lbl,'Position',[X0 Y 95 22]);
%         ui.(name) = uieditfield(setP,'numeric', ...
%             'Value',def,'Limits',[minv maxv], ...
%             'Position',[115 Y 90 22], ...
%             'ValueChangedFcn',@(src,~)fieldEdited(isScale));
%     end

    function addField(lbl,name,def,minv,maxv,isScale)
        if nargin<6, isScale = false; end
        Y = Y - dy;
        uilabel(setP,'Text',lbl,'Position',[X0 Y 95 22]);
        %         ui.(name) = uieditfield(setP,'numeric', ...
        %             'Value',def, ...
        %             'Limits',[minv maxv], ...          % ← now uses maxv you passed
        %             'ValueDisplayFormat','%.3g', ...   % nicer numbers
        %             'Position',[115 Y 90 22], ...
        %             'ValueChangedFcn',@(~,~)fieldEdited(isScale));
        
        ui.(name) = uieditfield(setP,'numeric', ...
            'Value',def,'Limits',[minv maxv], ...
            'ValueDisplayFormat','%.3g', ...   % 4 sig-figs
            'Position',[115 Y 90 22], ...
            'ValueChangedFcn',@(src,~)fieldEdited(isScale));
    end

    function fieldEdited(isScale)
        if isScale
            ui.autoScale = false;      % user takes control – stop autoscale
        end
        updatePlot();
    end

    function reflow
        W = fig.Position(3); H = fig.Position(4);
        ui.lst.Position = [margin margin+sliderH+margin listW H-3*margin-sliderH];
        setP.Position   = [W-panelW-margin margin+sliderH+margin panelW H-3*margin-sliderH];
        axX = margin + listW + margin;
        axW = W - axX - panelW - 2*margin;
        axH = H - 3*margin - sliderH;
        axY = margin + sliderH + margin;
        ui.ax.Position    = [axX axY axW axH];
        ui.sld.Position   = [axX margin axW sliderH-10];
        ui.curLabel.Position = [axX margin+sliderH-15 240 22];
    end
    function updatePlot()
        reflow();
        if ui.autoScale
            autoScaleNow();     % only while auto-mode is ON
        end
        t0=ui.timeEdit.Value; win=ui.lenEdit.Value;
        iBeg=max(1,floor(t0*fs)+1); iEnd=min(nsamples,floor((t0+win)*fs));
        chIdx=getCurrentIdx(); if isempty(chIdx),cla(ui.ax);return,end
        dat=ft_read_data(filename,'header',hdr,'begsample',iBeg,'endsample',iEnd,'chanindx',chIdx);
        if ui.fltChk.Value
            lo=ui.lowEdit.Value; hi=ui.highEdit.Value;
            if lo<hi && hi<fs/2, [b,a]=butter(4,[lo hi]/(fs/2)); dat=filtfilt(b,a,dat')'; end
        end
        
        for k=1:size(dat,1)
            lab=labNS{chIdx(k)}; scale=ui.othScaleEdit.Value;
            if endsWith(lab,'1'), scale=ui.magScaleEdit.Value;
            elseif endsWith(lab,{'2','3'}), scale=ui.gradScaleEdit.Value;
            elseif startsWith(lab,'EEG'), scale=ui.eegScaleEdit.Value;
            elseif startsWith(lab,'ECG'), scale=ui.ecgScaleEdit.Value; end
            dat(k,:)=dat(k,:)*scale;
        end
        sep=ui.gainEdit.Value*(std(dat(:),'omitnan')+eps);
        yOff=((size(dat,1)-1):-1:0)*sep; tvec=(iBeg:iEnd)/fs;
        cla(ui.ax);
        for k=1:size(dat,1), plot(ui.ax,tvec,dat(k,:)+yOff(k)); hold(ui.ax,'on'); end
        ui.ax.XLim=[t0 t0+win];
        if ui.classicChk.Value
            ui.ax.YTick=fliplr(yOff); ui.ax.YTickLabel=flip(labels(chIdx)); ui.ax.FontSize=10;
        else, ui.ax.YTick=[]; ui.ax.YTickLabel={}; end
        ui.curLabel.Text=sprintf('t = %.2f – %.2f s',t0,t0+win);
        
        if isgraphics(ui.selLine), delete(ui.selLine); end
        if ui.timeEdit.Value>=t0 && ui.timeEdit.Value<=t0+win
            ui.selLine=line(ui.ax,[ui.timeEdit.Value ui.timeEdit.Value],ui.ax.YLim,...
                'Color','k','LineWidth',1.5);
        end
        
        % -------- spike overlays ------------------------------------
        delete(findall(ui.ax,'Tag','SpikeLine'));
        spWin = ui.spikeT(ui.spikeT>=t0 & ui.spikeT<=t0+win);
        for tt = spWin(:).'
            line(ui.ax,[tt tt],ui.ax.YLim,'Color','m','LineStyle','--',...
                'Tag','SpikeLine','LineWidth',0.5);
        end
    end
    function autoScaleNow
        t0  = ui.timeEdit.Value;
        win = ui.lenEdit.Value;
        iBeg = max(1,floor(t0*fs)+1);
        iEnd = min(nsamples,floor((t0+win)*fs));
        
        chIdx = getCurrentIdx();
        if isempty(chIdx), return, end
        
        datRaw = ft_read_data(filename,'header',hdr, ...
            'begsample',iBeg,'endsample',iEnd, ...
            'chanindx',chIdx);
        chRMS  = std(datRaw,0,2,'omitnan');
        m      = @(sel) median(chRMS(sel),'omitnan');
        
        isMag  = endsWith(labNS(chIdx),'1');
        isGrad = endsWith(labNS(chIdx),'2') | endsWith(labNS(chIdx),'3');
        
        mMAG  = m(isMag);
        mGRAD = m(isGrad);
        mEEG  = m(startsWith(labNS(chIdx),'EEG'));
        mECG  = m(startsWith(labNS(chIdx),'ECG'));
        mOTH  = m(~(isMag | isGrad | startsWith(labNS(chIdx),{'EEG','ECG'})));
        
        medArray = [mMAG mGRAD mEEG mECG mOTH];
        hasVal   = ~isnan(medArray);
        
        if sum(hasVal) == 1
            ref = medArray(hasVal);            % • use that very RMS as reference
        else
            ref = median(medArray(hasVal));
        end
        % ------------------------------------------------------l------
        
        clip = @(v,lim) max(min(v,lim(2)),lim(1));
        
        digits = 3;                         % how many significant digits
        if mMAG  > 0
            ui.magScaleEdit.Value  = round( clip(ref/mMAG ,  ui.magScaleEdit.Limits), digits,'significant');
        end
        if mGRAD > 0
            ui.gradScaleEdit.Value = round( clip(ref/mGRAD, ui.gradScaleEdit.Limits), digits,'significant');
        end
        if mEEG  > 0
            ui.eegScaleEdit.Value  = round( clip(ref/mEEG , ui.eegScaleEdit.Limits), digits,'significant');
        end
        if mECG  > 0
            ui.ecgScaleEdit.Value  = round( clip(ref/mECG , ui.ecgScaleEdit.Limits), digits,'significant');
        end
        if mOTH  > 0
            ui.othScaleEdit.Value  = round( clip(ref/mOTH , ui.othScaleEdit.Limits), digits,'significant');
        end
        
    end
    function idx=getCurrentIdx()
        if ui.flvChk.Value, idx=1:numel(labels); ui.lst.Enable='off';
        elseif ~strcmp(ui.roiDrop.Value,'Manual')
            tgtNS=strrep(roiDefs(ui.roiDrop.Value),' ','');
            idx=find(ismember(labNS,tgtNS)); ui.lst.Enable='off';
        else, idx=find(ismember(labels,ui.lst.Value)); ui.lst.Enable='on';
        end
    end
    function openSpikePanel(ui, filename, hdr, fs, nSamples, getCurrentIdx, refreshMain)
        % openSpikePanel  Spike‑detection panel for the MNE Browser (compatible ≥R2018b)
        %
        %   This version detects whether the modern App‑Designer grid layout
        %   (uigridlayout, introduced in R2019b) is available.  If not, the panel
        %   falls back to absolute positioning so that it still works on older
        %   releases that support UIFigure but lack the newer layout managers.
        %
        %   ui            : struct from the main browser (stores spikeT etc.)
        %   filename/hdr  : FieldTrip raw file & header
        %   fs            : sampling rate of the original recording
        %   nSamples      : total # samples in the file
        %   getCurrentIdx : handle that returns the channel indices currently
        %                   selected in the browser
        %   refreshMain   : handle that forces a redraw of the main trace
        %
        %   VY 2025‑06‑30
        
        %%
        disp('spike being assessed ...')
        %% --------------------------------------------------------------------
        % Create the panel -----------------------------------------------------
        f = uifigure('Name','Spike Detection', ...
            'Position',[200 200 340 410], ...
            'Resize','off');
        
        % Close should simply hide/close the window and resume any uiwait calls
        f.CloseRequestFcn = @(src,~) delete(src);
        
        %% --------------------------------------------------------------------
        % Decide on layout method ---------------------------------------------
        useGrid = (exist('uigridlayout','file') == 2);   % true for R2019b+
        if useGrid
            gl            = uigridlayout(f,[7 3]);
            gl.RowHeight   = {22 22 22 22 32 '1x' 22};
            gl.ColumnWidth = {110 '1x' 90};
        end
        
        % Helper for absolute‑position fallback (pre‑R2019b) -------------------
        if ~useGrid
            colW = [110 140  90];           % widths of the 3 columns
            rowH = [22 22 22 22 32 200 22]; % heights of the 7 rows
            gap  = 4;                       % pixel gap between cells
            cellpos = @(r,c) [ ...          % anonymous that mimics grid cell
                10 + sum(colW(1:c-1)) + gap*(c-1), ...
                410 - sum(rowH(1:r)) - gap*r, ...
                colW(c), rowH(r) ];
            fullrow = @(r) [10, 410 - sum(rowH(1:r)) - gap*r, ...
                sum(colW) + gap*2, rowH(r) ];
        end
        
        %% --------------------------------------------------------------------
        % Row 1 – modality selector -------------------------------------------
        if useGrid
            uilabel(gl,'Text','Modality:','HorizontalAlignment','left');
            bg = uibuttongroup(gl,'Layout','row');
            uiradiobutton(bg,'Text','MEG','Value',true,'Tag','meg');
            uiradiobutton(bg,'Text','EEG','Tag','eeg');
            uilabel(gl,'Text','');   % spacer (3rd column)
        else
            uilabel(f,'Text','Modality:', 'HorizontalAlignment','left', ...
                'Position',cellpos(1,1));
            bg = uibuttongroup(f,'Position',cellpos(1,2));
            uiradiobutton(bg,'Text','MEG','Value',true,'Position',[5 2 60 18],'Tag','meg');
            uiradiobutton(bg,'Text','EEG','Position',[70 2 50 18],'Tag','eeg');
        end
        
        %% Row 2 – band‑pass ----------------------------------------------------
        if useGrid
            uilabel(gl,'Text','Band [Hz]:','HorizontalAlignment','left');
            bpField = uieditfield(gl,'text','Value','80 250');
            uilabel(gl,'Text','');   % spacer
        else
            uilabel(f,'Text','Band [Hz]:','HorizontalAlignment','left', ...
                'Position',cellpos(2,1));
            bpField = uieditfield(f,'text','Value','80 250', ...
                'Position',cellpos(2,2));
        end
        
        %% Row 3 – z‑threshold --------------------------------------------------
        if useGrid
            uilabel(gl,'Text','z‑score thr:','HorizontalAlignment','left');
            zField  = uieditfield(gl,'numeric','Value',5,'Limits',[1 50]);
            uilabel(gl,'Text','');
        else
            uilabel(f,'Text','z‑score thr:','HorizontalAlignment','left', ...
                'Position',cellpos(3,1));
            zField  = uieditfield(f,'numeric','Value',5,'Limits',[1 50], ...
                'Position',cellpos(3,2));
        end
        
        %% Row 4 – scope dropdown ----------------------------------------------
        if useGrid
            uilabel(gl,'Text','Scope:','HorizontalAlignment','left');
            scopeDrop = uidropdown(gl,'Items',{'Current window','Whole file'}, ...
                'Value','Current window');
            uilabel(gl,'Text','');
        else
            uilabel(f,'Text','Scope:','HorizontalAlignment','left', ...
                'Position',cellpos(4,1));
            scopeDrop = uidropdown(f,'Items',{'Current window','Whole file'}, ...
                'Value','Current window', ...
                'Position',cellpos(4,2));
        end
        
        %% Row 5 – RUN button ---------------------------------------------------
        if useGrid
            uilabel(gl,'Text','');              % empty first column (spacer)
            runBtn = uibutton(gl,'Text','Run','ButtonPushedFcn',@runDetect);
            uilabel(gl,'Text','');
        else
            runBtn = uibutton(f,'Text','Run','ButtonPushedFcn',@runDetect, ...
                'Position',cellpos(5,2));
        end
        
        %% Row 6 – spike table --------------------------------------------------
        if useGrid
            spkTbl = uitable(gl,'ColumnName',{'Time (s)'},'ColumnEditable',false);
        else
            spkTbl = uitable(f,'ColumnName',{'Time (s)'},'ColumnEditable',false, ...
                'Position',fullrow(6));
        end
        
        %% Row 7 – info label ---------------------------------------------------
        if useGrid
            infoLbl = uilabel(gl,'Text','','FontAngle','italic');
        else
            infoLbl = uilabel(f,'Text','', 'FontAngle','italic', ...
                'Position',fullrow(7));
        end
        
        %% --------------------------------------------------------------------
        % CACHED DATA -----------------------------------------------------------
        cache.raw   = [];   % resampled data segment
        cache.t0    = [];   % segment start (s)
        cache.fsRes = 200;  % resample Fs
        
        
        function runDetect(~,~)
            % =================================================================
            % 0) Progress-bar infrastructure
            % =================================================================
            useUIP = (exist('uiprogressdlg','file') == 2);   % >= R2019a
            dlg    = makeDlg(useUIP);                       % helper just below
            % close the dialog no matter how we exit
            dlgCleaner = onCleanup(@() safeDelete(dlg));
            
            % two small handles we use everywhere
            if useUIP
                setPct      = @(p,varargin) updateDlg(p,varargin{:});
                wasCanceled = @() dlg.CancelRequested;
            else
                setPct      = @(p,varargin) updateBar(p,varargin{:});
                wasCanceled = @() ~isvalid(dlg);            % waitbar closed
            end
            
            % =================================================================
            % 1) Parse UI
            % =================================================================
            setPct(0.00,'Parsing settings …');
            bp = sort(str2double(split(bpField.Value))');
            validateattributes(bp,{'numeric'}, ...
                {'numel',2,'finite','increasing','<=',fs/2});
            zthr  = zField.Value;
            whole = strcmp(scopeDrop.Value,'Whole file');
            modal = bg.SelectedObject.Tag;   % 'meg' or 'eeg'
            if wasCanceled(), return, end
            
            % =================================================================
            % 2) Obtain / cache raw data
            % =================================================================
            setPct(0.10,'Reading data …');
            if whole
                if isempty(cache.raw)
                    cache.raw = readRaw(1,nSamples);
                    cache.t0  = 0;
                end
                dat = cache.raw;   t0 = 0;
            else
                t0  = ui.timeEdit.Value;       % window start (s)
                win = ui.lenEdit.Value;
                s1  = max(1, floor(t0*fs)+1);
                s2  = min(nSamples, floor((t0+win)*fs));
                if isempty(cache.raw) || abs(cache.t0 - t0) > 1e-3
                    cache.raw = readRaw(s1,s2);
                    cache.t0  = t0;
                end
                dat = cache.raw;
            end
            if wasCanceled(), return, end
            
            % =================================================================
            % 3) Minimal preprocessing
            % =================================================================
            setPct(0.35,'Filtering …');
            if strcmp(modal,'meg')
                dat = detrend(single(dat)','constant')';
            end
            if bp(2) < fs/2
                d = designfilt('bandpassfir', ...
                    'StopbandFrequency1',bp(1)*0.8, ...
                    'PassbandFrequency1',bp(1), ...
                    'PassbandFrequency2',bp(2), ...
                    'StopbandFrequency2',bp(2)*1.2, ...
                    'SampleRate',fs);
                dat = filtfilt(d,double(dat'))';
            end
            if wasCanceled(), return, end
            
            % =================================================================
            % 4) (Optional) resample for speed
            % =================================================================
            setPct(0.55,'Resampling …');
            if fs ~= cache.fsRes
                dat  = resample(dat',cache.fsRes,fs)';   % channels × time
                rsFs = cache.fsRes;
            else
                rsFs = fs;
            end
            if wasCanceled(), return, end
            
            % =================================================================
            % 5) Spike detection (RMS z-score)
            % =================================================================
            setPct(0.75,'Detecting spikes …');
            env    = sqrt(movmean(dat.^2,0.03*rsFs,2));     % 30 ms RMS
            zenv   = (env - mean(env,2))./std(env,0,2);
            spikes = find(any(zenv > zthr,1));              % OR across chans
            spkT   = (unique(spikes)-1)/rsFs + t0;          % seconds (file)
            if wasCanceled(), return, end
            
            % =================================================================
            % 6) Update UI and main browser
            % =================================================================
            setPct(0.90,'Updating display …');
            ui.spikeT   = spkT;
            spkTbl.Data = num2cell(spkT(:));                % table update
            infoLbl.Text= sprintf('%d spike(s) detected',numel(spkT));
            refreshMain();
            if wasCanceled(), return, end
            
            % =================================================================
            % 7) Done
            % =================================================================
            setPct(1.00,'Done ✔');  pause(0.3);             % show 100 % briefly
            
            function h = makeDlg(useUIP)
                % Desired size for the progress window
                sz = [260 80];                               % [width height]
                
                if useUIP
                    % Modern progress dialog -----------------------------------------
                    h = uiprogressdlg(f, ...
                        'Title'    ,'Spike detection', ...
                        'Message'  ,'Initializing …', ...
                        'Cancelable','on');
                    
                    % Resize / center only if that MATLAB release supports it
                    if isprop(h,'Position')
                        parPos = f.Position;                                 % [x y w h]
                        newPos = [ ...
                            parPos(1) + (parPos(3)-sz(1))/2 , ...            % X-center
                            parPos(2) + (parPos(4)-sz(2))/2 , ...            % Y-center
                            sz ];
                        h.Position = newPos;
                    end
                    
                else
                    % Legacy waitbar --------------------------------------------------
                    h = waitbar(0,'Initializing …','Name','Spike detection', ...
                        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                    
                    parPos = f.Position;
                    newPos = [ ...
                        parPos(1) + (parPos(3)-sz(1))/2 , ...
                        parPos(2) + (parPos(4)-sz(2))/2 , ...
                        sz ];
                    set(h,'Units','pixels','Position',newPos);
                end
            end
            
            function updateDlg(p,msg)
                if nargin<2, msg=''; end
                dlg.Value   = p;
                dlg.Message = msg;
                drawnow;
            end
            
            function updateBar(p,msg)
                if nargin<2, msg=''; end
                waitbar(p,dlg,msg);
            end
            
            function safeDelete(h)
                if isvalid(h)
                    delete(h);
                end
            end
        end   % <–– END of runDetect
        
        % show 100 % briefly
        %         end
        
        % CALLBACKS -------------------------------------------------------------
        %         function runDetect(~,~)
        %             %–– parse UI ----------------------------------------------------
        %             bp = str2double(split(bpField.Value));
        %             bp = sort(bp(:)');
        %             validateattributes(bp,{'numeric'},{'numel',2,'finite','increasing','<=',fs/2});
        %             zthr  = zField.Value;
        %             whole = strcmp(scopeDrop.Value,'Whole file');
        %             modal = bg.SelectedObject.Tag;   % 'meg' or 'eeg'
        %
        %             %–– obtain / cache raw data ------------------------------------
        %             if whole
        %                 if isempty(cache.raw)
        %                     cache.raw = readRaw(1,nSamples);
        %                     cache.t0  = 0;
        %                 end
        %                 dat = cache.raw; t0 = 0;
        %             else
        %                 t0  = ui.timeEdit.Value;    % seconds at start of current window
        %                 win = ui.lenEdit.Value;
        %                 s1 = max(1, floor(t0*fs)+1);
        %                 s2 = min(nSamples, floor((t0+win)*fs));
        %                 if isempty(cache.raw) || abs(cache.t0 - t0) > 1e-3
        %                     cache.raw = readRaw(s1,s2);
        %                     cache.t0  = t0;
        %                 end
        %                 dat = cache.raw;
        %             end
        %
        %             %–– minimal preprocessing -------------------------------------
        %             if strcmp(modal,'meg')
        %                 dat = detrend(single(dat)','constant')';
        %             end
        %             if bp(2) < fs/2
        %                 d = designfilt('bandpassfir', ...
        %                     'StopbandFrequency1',bp(1)*0.8, ...
        %                     'PassbandFrequency1',bp(1), ...
        %                     'PassbandFrequency2',bp(2), ...
        %                     'StopbandFrequency2',bp(2)*1.2, ...
        %                     'SampleRate',fs);
        %                 dat = filtfilt(d,double(dat'))';
        %             end
        %
        %             %–– resample to lower rate (speed) -----------------------------
        %             if fs ~= cache.fsRes
        %                 dat  = resample(dat',cache.fsRes,fs)';
        %                 rsFs = cache.fsRes;
        %             else
        %                 rsFs = fs;
        %             end
        %
        %             %–– spike detection (RMS z‑score) ------------------------------
        %             env   = sqrt(movmean(dat.^2,0.03*rsFs,2));   % 30 ms RMS
        %             zenv  = (env - mean(env,2))./std(env,0,2);
        %             spikes = find(any(zenv > zthr,1));           % OR across channels
        %             spikes = unique(spikes);
        %             spkT   = (spikes-1)/rsFs + t0;               % seconds from file start
        %
        %             %–– update UI and main browser --------------------------------
        %             ui.spikeT = spkT;
        %             spkTbl.Data = num2cell(spkT(:));
        %             infoLbl.Text = sprintf('%d spike(s) detected',numel(spkT));
        %             refreshMain();
        %
        %         end
        
        function dat = readRaw(begSamp,endSamp)
            chIdx = getCurrentIdx();
            dat   = ft_read_data(filename,'header',hdr, ...
                'begsample',begSamp,'endsample',endSamp, ...
                'chanindx',chIdx);
        end
        
        disp('completed.')
        
    end

end
