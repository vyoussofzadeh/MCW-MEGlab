%% READING TRINGGER VALUES
function [detResp,detTrig] = do_readtriggers(cfg, datafile)

hdr = ft_read_header(datafile); %read header information
Fsample = hdr.Fs;


Index = strfind(hdr.label,{'STI101'});
Index = find(not(cellfun('isempty',Index)));

if isempty(Index)
    error('stimulus data, STI101, is missing from the data header, ');
end

detTrig=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','neuromag_fif','dataformat','neuromag_fif');
detTrig = (detTrig - min(detTrig));
detTrig=bitand(detTrig,255);

tt = linspace(1, length(detTrig)/Fsample, length(detTrig));

if cfg.plot ==1
    figure,
    subplot 211
    plot(tt,detTrig)
    title('detTrig')
end

detResp = ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','digital trigger','dataformat','neuromag_fif');
detResp = detResp - detTrig;
detResp = (detResp - min(detResp));
resp = (detResp - mean(detResp))./max(detTrig);

if cfg.plot ==1
    %     figure,
    subplot 212
    plot(tt,detResp)
    title('detResp')
    
    figure,
    plot(tt,detTrig)
    hold on
    plot(tt,resp, 'r')
    title('All triggers')
end

unique(detResp)
unique(detTrig)

end