function mne_browser_gui_V1_55(filename)
% MNE_BROWSER_GUI – Raw MEG/EEG browser (MATLAB UIFigure)
%
% v1.55 | 2025-06-15
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
        
        %         if mMAG  >0, ui.magScaleEdit.Value  = clip(ref/mMAG ,  ui.magScaleEdit.Limits); end
        %         if mGRAD >0, ui.gradScaleEdit.Value = clip(ref/mGRAD, ui.gradScaleEdit.Limits); end
        %         if mEEG  >0, ui.eegScaleEdit.Value  = clip(ref/mEEG , ui.eegScaleEdit.Limits); end
        %         if mECG  >0, ui.ecgScaleEdit.Value  = clip(ref/mECG , ui.ecgScaleEdit.Limits); end
        %         if mOTH  >0, ui.othScaleEdit.Value  = clip(ref/mOTH , ui.othScaleEdit.Limits); end
        
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
%     function autoScaleCb(~,~)
%         ui.autoScale = true;   % keep behaviour
%         autoScaleNow;
%         updatePlot;
%         t0=ui.timeEdit.Value; win=ui.lenEdit.Value;
%         iBeg=max(1,floor(t0*fs)+1); iEnd=min(nsamples,floor((t0+win)*fs));
%         chIdx=getCurrentIdx(); if isempty(chIdx),return,end
%         raw = ft_read_data(filename,'header',hdr,'begsample',iBeg,'endsample',iEnd,'chanindx',chIdx);
%         rms = std(raw,0,2,'omitnan');
%         med=@(sel)median(rms(sel),'omitnan');
%         isMag=endsWith(labNS(chIdx),'1'); isGrad=endsWith(labNS(chIdx),{'2','3'});
%         m=[med(isMag) med(isGrad) med(startsWith(labNS(chIdx),'EEG')) ...
%             med(startsWith(labNS(chIdx),'ECG')) med(~(isMag|isGrad|startsWith(labNS(chIdx),{'EEG','ECG'})))];
%         ref=median(m,'omitnan'); if isnan(ref)||ref==0,ref=1;end
%         clip=@(v,lim)max(min(v,lim(2)),lim(1));
%         if med(isMag)>0,  ui.magScaleEdit.Value = clip(ref/med(isMag),  ui.magScaleEdit.Limits); end
%         if med(isGrad)>0, ui.gradScaleEdit.Value= clip(ref/med(isGrad), ui.gradScaleEdit.Limits); end
%         if med(startsWith(labNS(chIdx),'EEG'))>0
%             ui.eegScaleEdit.Value = clip(ref/med(startsWith(labNS(chIdx),'EEG')), ui.eegScaleEdit.Limits); end
%         if med(startsWith(labNS(chIdx),'ECG'))>0
%             ui.ecgScaleEdit.Value = clip(ref/med(startsWith(labNS(chIdx),'ECG')), ui.ecgScaleEdit.Limits); end
%         if med(~(isMag|isGrad|startsWith(labNS(chIdx),{'EEG','ECG'})))>0
%             ui.othScaleEdit.Value = clip(ref/med(~(isMag|isGrad|startsWith(labNS(chIdx),{'EEG','ECG'}))), ui.othScaleEdit.Limits); end
%         updatePlot();
%     end
    function idx=getCurrentIdx()
        if ui.flvChk.Value, idx=1:numel(labels); ui.lst.Enable='off';
        elseif ~strcmp(ui.roiDrop.Value,'Manual')
            tgtNS=strrep(roiDefs(ui.roiDrop.Value),' ','');
            idx=find(ismember(labNS,tgtNS)); ui.lst.Enable='off';
        else, idx=find(ismember(labels,ui.lst.Value)); ui.lst.Enable='on';
        end
    end
    function openSpikePanel(ui,filename,hdr,fs,nsamples,getCurrentIdx,updatePlot)
        
        % ---------- MODAL DIALOG (classic figure, not uifigure) -------------
        d = dialog('Name','Spike Detection', ...
            'Position',[200 200 330 390], ...
            'WindowStyle','modal');
        
        % keep a struct where we remember user choices
        par = struct('modsel',1);   % 1 = MEG, 2 = EEG  (default MEG)
        
        % -------------------------- LAYOUT ----------------------------------
        % we build a tiny grid by hand (uigridlayout is not available here)
        ypos = 350; dy = 30; lblW = 110; fldW =  80;
        
        % ──-(A)  Modality selector ───────────────────────────────────────────
        uicontrol(d,'Style','text','String','Modality:', ...
            'Position',[10 ypos lblW 18], ...
            'HorizontalAlignment','left');
        
        % 1⃣  create the button-group first
        bg = uibuttongroup(d, ...
            'Position',[120 ypos-5 120 22], ...
            'SelectionChangedFcn',@selMod);   % ← nested cb, see below
        
        % 2⃣  the two radio-buttons (note the Parent property)
        uicontrol('Parent',bg,'Style','radiobutton','String','MEG', ...
            'Position',[0   0 50 18],'Tag','1','Value',1);
        uicontrol('Parent',bg,'Style','radiobutton','String','EEG', ...
            'Position',[60  0 50 18],'Tag','2');
        
        % ── nested callback (must sit **inside** openSpikePanel) ─────────────
        function selMod(~,ev)
            par.modsel = str2double(ev.NewValue.Tag);    % 1 = MEG, 2 = EEG
        end
        
        % (B) –– Band-pass field -------------------------------------------
        ypos = ypos-dy;
        uicontrol(d,'Style','text','String','Band [Hz]:', ...
            'Position',[10 ypos lblW 18],'HorizontalAlignment','left');
        bpField = uicontrol(d,'Style','edit','String','4 40', ...
            'Position',[120 ypos fldW 22]);
        
        % (C) –– z-threshold -------------------------------------------------
        ypos = ypos-dy;
        uicontrol(d,'Style','text','String','z-score thr:', ...
            'Position',[10 ypos lblW 18],'HorizontalAlignment','left');
        zField  = uicontrol(d,'Style','edit','String','7', ...
            'Position',[120 ypos fldW 22]);
        
        % (D) –– scope (current / whole) ------------------------------------
        ypos = ypos-dy;
        uicontrol(d,'Style','text','String','Scope:', ...
            'Position',[10 ypos lblW 18],'HorizontalAlignment','left');
        scope = uicontrol(d,'Style','popupmenu', ...
            'String',{'Current window','Whole file'}, ...
            'Position',[120 ypos fldW+40 22]);
        
        % (E) –– RUN button --------------------------------------------------
        uicontrol(d,'Style','pushbutton','String','Run', ...
            'Position',[120 ypos-40 60 25], ...
            'Callback',@doDetect);
        
        % (F) –– spike list table -------------------------------------------
        uicontrol(d,'Style','text','String','Spikes (s):', ...
            'FontWeight','bold','HorizontalAlignment','left', ...
            'Position',[10 100 100 18]);
        spikeTbl = uitable(d,'Position',[10 10 300 90], ...
            'ColumnName',{'Time (s)'}, ...
            'CellSelectionCallback',@spikeSelect);
        
        % --------------------  DETECT CALLBACK  ----------------------------
        function doDetect(~,~)
            bp    = sort(str2num(bpField.String));   %#ok<ST2NM>
            zthr  = str2double(zField.String);
            whole = scope.Value==2;
            
            % read data
            if whole
                dat = ft_read_data(filename,'header',hdr);
            else
                % current view window only
                t0   = ui.timeEdit.Value;
                win  = ui.lenEdit.Value;
                iBeg = max(1,floor(t0*fs)+1);
                iEnd = min(nsamples,floor((t0+win)*fs));
                chIdx= getCurrentIdx();
                dat  = ft_read_data(filename,'header',hdr, ...
                    'begsample', iBeg, 'endsample', iEnd, ...
                    'chanindx',  chIdx);
            end
            
            % 1. --------- read & (optionally) preprocess --------------------
            switch par.modsel
                case 1      % ---------- MEG ----------
                    cfg = []; cfg.dataset = filename;
                    cfg.channel = {'meg','eog'};         % mags+grads+EOG
                    raw_data = ft_preprocessing(cfg);
                    modal    = 'meg*';   cutoff_val = 10;
                otherwise   % ---------- EEG ----------
                    cfg = []; cfg.dataset = filename;
                    cfg.channel = {'eeg','eog'};
                    cfg.reref   = 'yes'; cfg.refmethod='avg'; cfg.refchannel='all';
                    raw_data = ft_preprocessing(cfg);
                    modal    = 'eeg*';   cutoff_val = 5;
            end
            
            % optional resample (kept from your snippet)
            cfg = []; cfg.resamplefs = 200;
            dataForSpk = ft_resampledata(cfg, raw_data);
            
            % 2. --------- build cfg for detMEGspikes ------------------------
            cfgdet = struct('dtype', modal(1:3), ...
                'fs',    dataForSpk.fsample, ...
                'spikeband',bp, ...
                'zthresh',zthr, ...
                'minterval',0.03, ...
                'statswin',2, ...
                'revspikes',0, ...
                'statswinType','global');
            
            [~, allSpk] = detMEGspikes(dataForSpk.trial{1}, string(dataForSpk.label), cfgdet);
            ui.spikeT   = allSpk ./ dataForSpk.fsample;   % → seconds
            
            % 3. --------- update table & main plot -------------------------
            spikeTbl.Data = num2cell(ui.spikeT(:));
            updatePlot();
            msgbox(sprintf('%d spikes detected',numel(ui.spikeT)),'Spike detection');
        end
    end

%     function openSpikePanel(ui, filename, hdr, fs, nSamples, getCurrentIdx, refreshMain)
%         % Optimised Spike-detection panel for MNE-Browser 1.55+
%         %
%         % ui : handle to browser's UI struct (gets spikeT, cache fields)
%         % refreshMain : function handle to redraw the main trace
%         
%         %% ----------  FIGURE ----------------------------------------------------
%         % f = uifigure( ...
%         %         'Name','Spike Detection', ...
%         %         'Position',[200 200 340 410], ...
%         %         'WindowStyle','modal', ...           % keeps focus but main UI alive
%         %         'Resize','off');
%         
%         f = uifigure('Name','Spike Detection', ...
%             'Position',[200 200 340 410], ...
%             'Resize','off');
%         
%         movegui(f,'center')        % nice-to-have
%         
%         % block MATLAB and the parent UI until the window closes
%         uiwait(f);                  % ↓ in the CloseRequestFcn call uiresume(f)
%         
%         f.CloseRequestFcn = @(src,~) uiresume(src);
%         
%         gl = uigridlayout(f,[7 3]);
%         gl.RowHeight   = {22 22 22 22 32 '1x' 'fit'};
%         gl.ColumnWidth = {110 '1x' 90};
%         
%         %% ----------  CONTROLS --------------------------------------------------
%         % Row 1 – modality
%         % uilabel(gl,'Text','Modality:',                 'HorizontalAlignment','left');
%         % bg = uibuttongroup(gl,'Layout','row');
%         % uiradiobutton(bg,'Text','MEG','Value',true,'Tag','meg');
%         % uiradiobutton(bg,'Text','EEG','Tag','eeg');
%         % % uigridcontainer(gl);        % empty placeholder
%         % uilabel(gl,'Text','');
%         
%         uilabel(gl,'Text','Modality:','HorizontalAlignment','left');
%         bg = uibuttongroup(gl,'Layout','row');
%         uiradiobutton(bg,'Text','MEG','Value',true,'Tag','meg');
%         uiradiobutton(bg,'Text','EEG','Tag','eeg');
%         uilabel(gl,'Text','');   % ← replaces uigridcontainer
%         
%         % Row 2 – band-pass
%         uilabel(gl,'Text','Band [Hz]:','HorizontalAlignment','left');
%         bpField = uieditfield(gl,'text','Value','80 250');
%         uigridcontainer(gl);
%         
%         % Row 3 – z-threshold
%         uilabel(gl,'Text','z-score thr:','HorizontalAlignment','left');
%         zField  = uieditfield(gl,'numeric','Value',5,'Limits',[1 50]);
%         uigridcontainer(gl);
%         
%         % Row 4 – scope
%         uilabel(gl,'Text','Scope:','HorizontalAlignment','left');
%         scopeDrop = uidropdown(gl, ...
%             'Items',{'Current window','Whole file'},'Value','Current window');
%         uigridcontainer(gl);
%         
%         % Row 5 – RUN button
%         runBtn = uibutton(gl,'Text','Run','ButtonPushedFcn',@runDetect);
%         
%         % Row 6 – spike table
%         spkTbl = uitable(gl,'ColumnName',{'Time (s)'},'ColumnEditable',false);
%         
%         % Row 7 – info label
%         infoLbl = uilabel(gl,'Text','','FontAngle','italic');
%         
%         %% ----------  CACHED DATA ----------------------------------------------
%         cache.raw   = [];   % resampled data segment
%         cache.t0    = [];   % segment start (s)
%         cache.fsRes = 200;  % resample Fs
%         
%         %% ----------  CALLBACKS -------------------------------------------------
%         function runDetect(~,~)
%             %–––– parse UI
%             bp    = str2double(split(bpField.Value));   bp = sort(bp(:)');
%             validateattributes(bp,{'numeric'},{'numel',2,'finite','increasing','<=',fs/2});
%             zthr  = zField.Value;
%             whole = strcmp(scopeDrop.Value,'Whole file');
%             modal = bg.SelectedObject.Tag;            % 'meg' / 'eeg'
%             
%             %–––– get (or reuse) raw chunk
%             if whole
%                 if isempty(cache.raw)
%                     cache.raw = readRaw(1, nSamples);
%                     cache.t0  = 0;
%                 end
%                 dat = cache.raw;
%                 t0  = 0;
%             else
%                 t0  = ui.timeEdit.Value;
%                 win = ui.lenEdit.Value;
%                 s1  = max(1 , floor(t0*fs)+1);
%                 s2  = min(nSamples, floor((t0+win)*fs));
%                 if isempty(cache.raw) || abs(cache.t0 - t0) > 1e-3
%                     cache.raw = readRaw(s1, s2);
%                     cache.t0  = t0;
%                 end
%                 dat = cache.raw;
%             end
%             
%             %–––– minimal preprocessing
%             if modal == "meg"
%                 dat = detrend(single(dat)','constant')';     % remove DC
%             end
%             if bp(2) < fs/2
%                 % Parks-McClellan linear phase FIR (order ≈ fs)
%                 d = designfilt('bandpassfir','StopbandFrequency1',bp(1)*0.8, ...
%                     'PassbandFrequency1',bp(1), ...
%                     'PassbandFrequency2',bp(2), ...
%                     'StopbandFrequency2',bp(2)*1.2, ...
%                     'SampleRate',fs);
%                 dat = filtfilt(d,double(dat') )';
%             end
%             %–––– resample
%             if fs ~= cache.fsRes
%                 dat = resample(dat', cache.fsRes, fs)';
%                 rsFs = cache.fsRes;
%             else
%                 rsFs = fs;
%             end
%             
%             %–––– spike detection (simple RMS-zscore)
%             env   = sqrt(movmean(dat.^2, 0.03*rsFs, 2));   % 30-ms RMS
%             zenv  = (env - mean(env,2))./std(env,0,2);
%             spikes = find(any(zenv > zthr,1));             % logical OR across channels
%             spikes = unique(spikes);                       % remove duplicates
%             spkT   = (spikes-1)/rsFs + t0;                 % seconds in file
%             
%             %% update UI + browser
%             ui.spikeT = spkT;
%             spkTbl.Data = num2cell(spkT(:));
%             infoLbl.Text = sprintf('%d spikes detected', numel(spkT));
%             refreshMain();
%         end
%         
%         function dat = readRaw(begSamp, endSamp)
%             chIdx = getCurrentIdx();
%             dat   = ft_read_data(filename,'header',hdr, ...
%                 'begsample',begSamp,'endsample',endSamp, ...
%                 'chanindx',chIdx);
%         end
%     end


end
