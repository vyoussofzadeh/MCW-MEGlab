function    [ttrl, trl1, trl_Info, evt_short, detTrig, resp]  = trialfun_sentence_meaningfulness_SP(cfg)

% trial function for the sentence meaningfulness task
% Author: Vahab Youssof Zadeh, vyoussofzadeh@mcw.edu
% Update date: 08/27/2021

%%
plotflag = cfg.plotflag;
datafile = cfg. datafile;
prestimTime = cfg.trialdef.prestimTime;
poststimTime = cfg.trialdef.poststimTime;

switch cfg.oneset_type
    case 'sentence'
        onsval = 5;
    case 'fixation'
        onsval = 6;
    case 'target'
        onsval = 7;
    case 'response'
        onsval = 8;
end

%%
hdr = ft_read_header(datafile); %read header information
Fsample = hdr.Fs;
prestimSamples = floor(prestimTime*Fsample);
poststimSamples = floor(poststimTime*Fsample);

%%
Index = strfind(hdr.label,{'STI101'});
Index = find(not(cellfun('isempty',Index)));

if isempty(Index)
    error('stimulus data, STI101, is missing from the data header, ');
end

%%
detTrig=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','neuromag_fif','dataformat','neuromag_fif');
detTrig = (detTrig - min(detTrig));
detTrig=bitand(detTrig,255);

if plotflag ==1
    figure,
    subplot 211
    plot(detTrig)
    title('detTrig')
end

detResp=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','digital trigger','dataformat','neuromag_fif');
detResp = detResp - detTrig;
detResp = (detResp - min(detResp));
resp = (detResp - mean(detResp))./max(detTrig);

if plotflag ==1
    %     figure,
    subplot 212
    plot(detResp)
    title('detResp')
end

%%
psample = 20; % compensation samples due to chnages in trigger rising and falling edges

%% Trigger noise reduction
idx = find(abs(diff(detTrig))  > 10);
detTrig1 = detTrig;
for i=1:length(idx)
    detTrig1(idx(i)-psample:idx(i)) = detTrig1(idx(i)-psample);
    detTrig1(idx(i):idx(i)+psample) = detTrig1(idx(i)+psample);
end
detTrig = detTrig1;

%%
if plotflag ==1
    figure,
    plot(detTrig)
    hold on
    plot(resp, 'r')
    title('All triggers')
end

%% Identify triggers, Up (rising  edge)-Trig, Down (falling  edge)-trig and Up (rising edge)-responses,

ixp_pos = find(diff(detTrig) > 0 ) + 1;
trgp_pos = detTrig(ixp_pos);

neg_ixp = find(diff((detTrig)) < 0 ) + 1;
neg_trgp = detTrig(neg_ixp);

ixpr = find(diff(abs(resp)) > 0 ) + 1;
trgpr = detTrig(ixpr);

%% Create trlBlockInfo
neweventsample = neg_ixp; eventval = neg_trgp;
cfg1 = [];
cfg1.detTrig = detTrig;
cfg1.neweventsample = neweventsample;
cfg1.eventval = eventval;
cfg1.plotflag = plotflag;
cfg1.psample = psample;
trl_neg = trig_detect1(cfg1);

neweventsample = ixp_pos; eventval = trgp_pos;
cfg1.detTrig = detTrig;
cfg1.neweventsample = neweventsample;
cfg1.eventval = eventval;
cfg1.plotflag = plotflag;
cfg1.psample = psample;
trl_pos = trig_detect1(cfg1);

% figure,plot(abs(resp))
neweventsample = ixpr; eventval = trgpr;
cfg1.resp = resp;
cfg1.dresp = diff(resp);
cfg1.neweventsample = neweventsample;
cfg1.eventval = eventval;
cfg1.plotflag = plotflag;
cfg1.psample = psample;
trl_resp = trig_resp(cfg1);

%% %All EVents: combining (and sorting) pos, neg, and response triggers (based on the sample number)
trl_all = [trl_resp; trl_pos ;trl_neg];
trl_all(isnan(trl_all(:,2))==1,:)=[];
[~, idx] = sort(trl_all(:,1));
trl_all_sort = trl_all(idx,:);

%% Remove redundant (sequential) sample triggers (mixed-up edges).
% df = diff(trl_all_sort(:,1));
trl = trl_all_sort;

%% Update/incorporate response accuracy, consider no-answer as correct
trl1 = trl;
A = sort(unique(trl1(:,5)))';
indices = strfind(diff(A), [1 1 1]);
if length(indices) > 1
    if A(indices(1)) > A(indices(2))
        indices_Sel = indices(1);
    else
        indices_Sel = indices(2);
    end
    indices = indices_Sel;
end

if ~isempty(indices)
    Cor_val = [A(indices),A(indices)+1, A(indices)+2, A(indices)+3];
    inCor_val = Cor_val-16;
    idx_targettype = ~isnan(trl1(:,3))==1;
    idx = find(idx_targettype ==1);
    for i=1:length(idx)
        if idx(i) > 2
            if ~isempty(find(trl1(idx(i),5) == inCor_val, 1)) ...
                    && ((trl1(idx(i)-1,2)==8) || (trl1(idx(i)-2,2)==8)) ...
                    && trl1(idx(i),2) < 10
                trl1(idx(i),2) = trl1(idx(i),2) + trl1(idx(i),2)*10;
            end
        end
    end
else
    indices = strfind(diff(A), [1]);
    if ~isempty(indices)
        inCor_val = [];
        Cor_val = [A(indices(1)),A(indices(1))+1, A(indices(2)),A(indices(2))+1];
        idx_targettype = ~isnan(trl1(:,3))==1;
        idx = find(idx_targettype ==1);
        for i=1:length(idx)
            if idx(i) > 2
                if ~isempty(find(trl1(idx(i),5) == inCor_val, 1)) ...
                        && ((trl1(idx(i)-1,2)==8) || (trl1(idx(i)-2,2)==8)) ...
                        && trl1(idx(i),2) < 10  && trl1(idx(i)-1,5) > 30
                    trl1(idx(i),2) = trl1(idx(i),2) + trl1(idx(i),2)*10;
                    %                     disp(i)
                end
            end
        end
    end
end

trl2 = trl1;
idx_targettype = ~isnan(trl2(:,3))==1;
idx = find(idx_targettype ==1);
for i=1:length(idx)
    if idx(i) > 2
        if ~isempty(find(trl2(idx(i),5) == inCor_val, 1)) ...
                && ((trl2(idx(i)-1,2)~=8) && (trl2(idx(i)-2,2)~=8)) ...
                && (trl2(idx(i)-1,2)==7) ...
                && trl2(idx(i),2) > 10
            %             disp(i)
            trl2(idx(i),2) = trl2(idx(i),2) - floor(trl2(idx(i),2)/10)*10;
        end
    end
end

%% evt_short (for quick testing triggers etc.)
idx_targetonset = trl2(:,2) == onsval;
idx_targetonset1 = trl2(:,5) == mode(trl2(idx_targetonset,5));
event_val = trl2(idx_targettype,2);
event_smpl = trl2(idx_targetonset1,1);
if length(event_smpl) == length(event_val)
    evt_short = [[1:length(event_smpl)]', event_smpl, event_val];
else
    trl3 = trl2;
    trl3(isnan(trl3(:,2))==1,:)=[];
    idx = (find(trl3(:,2) == onsval));
    idx2 = (find(~isnan(trl3(:,3))==1));
    a = zeros(60,2);
    a(1:length(idx),1) = idx;
    a(1:length(idx2),2) = idx2;
    dff = a(:,2) - a(:,1); idx1 = find(abs(dff) < 10);
    %     disp(a)
    a(idx1,2)= a(idx1,2)-dff(idx1);
    a  = a(1:min(length(idx), length(idx2)),:);
    [~,IA,IB] =  intersect(a(:,1),a(:,2));
    evt_short = [[1:length(IA)]', event_smpl(IA), event_val(IB)];
    warning('check the triggers, somthing is wrong!');
end

%%
ttrl = [];
ttrl(:,1) = evt_short(:,2) - prestimSamples;
ttrl(:,2) = evt_short(:,2) + poststimSamples;
ttrl(:,3) = - prestimSamples.*ones(size(evt_short,1),1);
ttrl(:,4) = evt_short(:,3);

%%
trl_Info={
    '1. Sample number'
    '2. Trigger value: 1. expected, 2. unexpected, 3. anomalous, 4. pseudoword, 5. sentence onset, 6, fixation onset, 7. target onset, 8. response onset'
    '3. Target occurrence, 1. yes, Nan. no'
    '4. Button press (response), left/right (1/2), '
    '5. Actual trigger value (from the header file after reading and bitand operation), '
    };

end

function trl = trig_detect1(cfg1)

neweventsample = cfg1.neweventsample;
eventval = cfg1.eventval;
plotflag = cfg1.plotflag;
psample = cfg1.psample;
detTrig = cfg1.detTrig;

trl = [];
for i=1:length(neweventsample)
    
    if plotflag ==1
        hold on
        y = ylim; % current y-axis limits
        plot([neweventsample(i) neweventsample(i)],y)
    end
    
    trl(i,1) =  neweventsample(i);
    eventval(i) = detTrig(neweventsample(i)+psample);
    
    switch eventval(i) % target type
        case {33, 41, 17, 25}
            trl(i,2) = 1; trl(i,3) = 1; % Expected - correct
        case {34, 42, 18, 26}
            trl(i,2) = 2; trl(i,3) = 1;  % Unexpected - correct
        case {35, 43, 19, 27}
            trl(i,2) = 3; trl(i,3) = 1;  % Anam - correct
        case {36, 44, 20, 28}
            trl(i,2) = 4; trl(i,3) = 1;  % Pseudo - correct
        case {9}
            trl(i,2) = 1; trl(i,3) = 1;  % Expected - incorrect
        case {10}
            trl(i,2) = 2; trl(i,3) = 1;  % Unexpected - incorrect
        case {11}
            trl(i,2) = 3; trl(i,3) = 1;  % Anam - incorrect
        case {12}
            trl(i,2) = 4; trl(i,3) = 1;  % Pseudo - incorrect
        case {64, 56}
            trl(i,2) = 5; trl(i,3) = nan;  % Sent onset - nan
        case {128, 120}
            trl(i,2) = 6; trl(i,3) = nan;  % Fixation onset - nan
        case {8, 0}
            trl(i,2) = 7; trl(i,3) = nan;  % Target onset - nan
        otherwise
            trl(i,2) = nan; trl(i,3) = nan;  % nan - nan
    end
    trl(i,4) = nan;
    trl(i,5) = eventval(i);
end
end

function trl = trig_resp(cfg1)

neweventsample = cfg1.neweventsample;
eventval = cfg1.eventval;
plotflag = cfg1.plotflag;
psample = cfg1.psample;
resp = cfg1.resp;
dres = cfg1.dresp;

% resp(abs(resp) < 0.5.*max(resp)) = 0;

trl = [];
for i=1:length(neweventsample)
    
    if plotflag ==1
        hold on
        y = ylim; % current y-axis limits
        plot([neweventsample(i) neweventsample(i)],y)
    end
    trl(i,1) =  neweventsample(i);
    
    eventval(i) = resp(neweventsample(i)-psample);
    
    A = resp(neweventsample(i)+psample);
    B = resp(neweventsample(i)-psample);
    if abs(A) > abs(B)
        eventval(i) = A;
    else
        eventval(i) = B;
    end
    respval_before = resp(neweventsample(i) -psample); % response values, trigger down
    respval_after = resp(neweventsample(i) + psample); % response values, triggerup
    
    %- responses, left and right
    if (i==1) && (abs(respval_before) > 50 || abs(respval_after) > 50)
        trl(i,2) = 8;
        if abs(respval_before) > abs(respval_after)
            trl(i,4) = 1; trl(i,3) = nan;
        else
            trl(i,4) = 2; trl(i,3) = nan;
        end
    else
        trl(i,2) = nan; trl(i,3) = nan; trl(i,4) = nan;
    end
    
    %- responses, left and right
    if (i>1) && (abs(respval_before) > 50 || abs(respval_after) > 50) && ((neweventsample(i) - neweventsample(i-1) > 500) || eventval(i) - eventval(i-1) > 30)
        trl(i,2) = 8;
        if abs(respval_before) > abs(respval_after)
            trl(i,4) = 1; trl(i,3) = nan;
        else
            trl(i,4) = 2; trl(i,3) = nan;
        end
    else
        trl(i,2) = nan; trl(i,3) = nan; trl(i,4) = nan;
    end
    trl(i,5) = eventval(i);
end
end