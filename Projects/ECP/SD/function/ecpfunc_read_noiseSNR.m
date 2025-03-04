function ecpfunc_read_noiseSNR()
% run_read_noiseSNR.m
%
% Example script for reading SNR data from multiple .mat files
% and incorporating them into your LI analysis or other pipelines.

    % -- (A) Identify all relevant .mat files in your directory -----------
    % Adjust the 'snr_*.mat' pattern if your filenames follow a different convention
    snrFiles = dir('snr_*.mat');
    
    % If you'd rather load each file explicitly, you can define them by hand:
    % snrFiles = {
    %     'snr_grad_tSSSvsRaw_1.mat'
    %     'snr_grad_MEGnet_Vs_tSSS_2.mat'
    %     'snr_mag_MEGnet_Vs_tSSS_2.mat'
    %     'snr_mag_tSSSvsRaw_1.mat'
    %     'snr_mag_tSSSvsRaw_2.mat'
    % };
    
    % -- (B) Initialize containers for the SNR data -----------------------
    all_snr_grad = [];  % e.g. you could store SNR grad data in a struct array
    all_snr_mag  = [];  % similarly for SNR mag
    
    % -- (C) Loop through each file and load data -------------------------
    for i = 1:length(snrFiles)
        fprintf('Loading file: %s\n', snrFiles(i).name);
        
        % Load the MAT file into a temporary struct
        tmp = load(snrFiles(i).name);
        
        % Check for presence of snr_grad or snr_mag structure
        if isfield(tmp, 'snr_grad')
            % If the file has a .snr_grad field
            snr_grad_data = tmp.snr_grad;
            
            % Example of how you might store each files data
            all_snr_grad(end+1).filename = snrFiles(i).name;  %#ok<AGROW>
            all_snr_grad(end).sub       = snr_grad_data.sub;
            all_snr_grad(end).sub_run   = snr_grad_data.sub_run;
            all_snr_grad(end).F_all     = snr_grad_data.F_all;
            all_snr_grad(end).snr_value = snr_grad_data.snr_value;

        elseif isfield(tmp, 'snr_mag')
            % If the file has a .snr_mag field
            snr_mag_data = tmp.snr_mag;
            
            all_snr_mag(end+1).filename = snrFiles(i).name;   %#ok<AGROW>
            all_snr_mag(end).sub       = snr_mag_data.sub;
            all_snr_mag(end).sub_run   = snr_mag_data.sub_run;
            all_snr_mag(end).F_all     = snr_mag_data.F_all;
            all_snr_mag(end).snr_value = snr_mag_data.snr_value;
            
        else
            warning('No recognizable SNR struct in %s (skipping).', snrFiles(i).name);
        end
    end
    
    % -- (D) Now incorporate the loaded SNR data into your LI analysis ----
    % Below is just a placeholder illustrating where you would integrate
    % whatever function or analysis you do for LI.

    % Example of calling an LI function for each subjects data:
    %
    %   for s = 1:length(all_snr_grad)
    %       current_snr = all_snr_grad(s).snr_value;
    %       current_sub = all_snr_grad(s).sub;
    %       % ... do your analysis with current_snr ...
    %
    %       % Example: lateralization_index = compute_LI(current_snr);
    %       % (Replace "compute_LI" with your real function!)
    %   end
    
    % If you have a single function to which you pass all data:
    %   results = compute_LI_all(all_snr_grad, all_snr_mag);
    
    % For demonstration, just show how many entries we loaded:
    fprintf('\nFinished loading:\n');
    fprintf('  SNR grad entries: %d\n', length(all_snr_grad));
    fprintf('  SNR mag entries:  %d\n', length(all_snr_mag));
    
    % -- (E) Save or return your aggregated data for further processing ---
    % If desired, you can save the aggregated data to a new .mat:
    % save('all_snr_data.mat', 'all_snr_grad', 'all_snr_mag', '-v7.3');
    
    % Or just end the function here; now you have loaded data in your workspace
    % if you run this as a script rather than a function.
end
