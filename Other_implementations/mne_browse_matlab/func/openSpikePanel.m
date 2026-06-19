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

% (E) –– time interval -------------------------------------------------
ypos = ypos - dy;
uicontrol(d,'Style','text','String','From (s):', ...   % <-- NEW
          'Position',[10 ypos lblW 18],'HorizontalAlignment','left');
tBegField = uicontrol(d,'Style','edit','String','', ...% <-- NEW
          'Position',[120 ypos fldW 22]);

ypos = ypos - dy;
uicontrol(d,'Style','text','String','To   (s):', ...
          'Position',[10 ypos lblW 18],'HorizontalAlignment','left');
tEndField = uicontrol(d,'Style','edit','String','', ...
          'Position',[120 ypos fldW 22]);



% % (To) field  ..........................................................
% ypos = ypos-dy;
% uicontrol(d,'Style','text','String','To      (s):', ...
%           'Position',[10 ypos lblW 18],'HorizontalAlignment','left');
% tEndField = uicontrol(d,'Style','edit','String','', ...
%           'Position',[120 ypos fldW 22]);
% 
% % >>> add this: drop by one row so we don’t paint over the boxes
% ypos = ypos - dy;
% 
% % (Run) button  ........................................................
% uicontrol(d,'Style','pushbutton','String','Run', ...
%           'Position',[120 ypos 60 25], ...
%           'Callback',@doDetect);

% (E) –– RUN button --------------------------------------------------
ypos = ypos - dy;
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
        
        % ------------------------------------------------ choose time window
if scope.Value == 2         % “Whole file”
    tBeg = 0;
    tEnd = nsamples/fs;
else                        % “Current / custom”
    % if user left boxes empty fall back to current view
    if isempty(tBegField.String)
        tBeg = ui.timeEdit.Value;
    else
        tBeg = str2double(tBegField.String);
    end
    if isempty(tEndField.String)
        tEnd = tBeg + ui.lenEdit.Value;
    else
        tEnd = str2double(tEndField.String);
    end
end

% sanity-clamp & convert to sample indices
tBeg = max(0,            tBeg);
tEnd = min(nsamples/fs,  tEnd);
if tEnd <= tBeg
    errordlg('"To" must be larger than "From"'); return
end
iBeg = floor(tBeg*fs) + 1;
iEnd = floor(tEnd*fs);

        
        
%         % ---------------------------------------------------------------------
% % 0) decide what time range we will read
% if whole                         % --- “Whole file” selected
%     iBeg = 1;
%     iEnd = nsamples;
% else                             % --- obey t-begin / t-end if given
%     % default to current window if boxes are empty
%     if isempty(tBegField.String), tBeg = ui.timeEdit.Value;        %#ok<*NODEF>
%     else,                         tBeg = str2double(tBegField.String); end
%     if isempty(tEndField.String), tEnd = tBeg + ui.lenEdit.Value;
%     else,                         tEnd = str2double(tEndField.String); end
% 
%     % safety clip
%     tBeg = max(tBeg ,0);
%     tEnd = min(tEnd , nsamples/fs);
%     if tEnd <= tBeg, errordlg('“To” must be > “From”'); return; end
% 
%     iBeg = floor(tBeg*fs)+1;
%     iEnd = floor(tEnd*fs);
% end

% finally read that chunk, **all visible channels**
chIdx = getCurrentIdx();
dat   = ft_read_data(filename,'header',hdr, ...
                     'begsample',iBeg,'endsample',iEnd, ...
                     'chanindx', chIdx);
% % ---------------------------------------------------------------------
% 
%         
%         % read data
%         if whole
%             dat = ft_read_data(filename,'header',hdr);
%         else
%             % current view window only
%             t0   = ui.timeEdit.Value;
%             win  = ui.lenEdit.Value;
%             iBeg = max(1,floor(t0*fs)+1);
%             iEnd = min(nsamples,floor((t0+win)*fs));
%             chIdx= getCurrentIdx();
%             dat  = ft_read_data(filename,'header',hdr, ...
%                 'begsample', iBeg, 'endsample', iEnd, ...
%                 'chanindx',  chIdx);
%         end
        
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
