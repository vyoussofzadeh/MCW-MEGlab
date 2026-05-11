function data = do_appyaft (cfg_main, data)

fs = data.fsample;
aft = cfg_main.aft;
aftval = cfg_main.aftval;

val = [];
val.time = data.time{1};
val.pow = data.trial{1};

for i=1:size(aft.eog,1)
    [~, idx1] = min(abs(val.time - aft.eog(i,1)/fs));
    [~, idx2] = min(abs(val.time - aft.eog(i,2)/fs));
    val.pow(:,idx1:idx2) = aftval;
end
for i=1:size(aft.jump,1)
    [~, idx1] = min(abs(val.time - aft.jump(i,1)/fs));
    [~, idx2] = min(abs(val.time - aft.jump(i,2)/fs));
    val.pow(:,idx1:idx2) = aftval;
end
for i=1:size(aft.rejseg,1)
    [~, idx1] = min(abs(val.time - aft.rejseg(i,1)/fs));
    [~, idx2] = min(abs(val.time - aft.rejseg(i,2)/fs));
    val.pow(:,idx1:idx2) = aftval;
end

data.time{1} = val.time;
data.trial{1} = val.pow;

end