function [all_snr_grad, all_snr_mag] = ecpfunc_read_noiseSNR()
% ecpfunc_read_noiseSNR
%
% Example script for reading SNR data from multiple .mat files
% and incorporating them into your LI analysis or other pipelines.
%
% USAGE:
%   [all_snr_grad, all_snr_mag] = ecpfunc_read_noiseSNR();
%
% DESCRIPTION:
%   This function loads any .mat file matching the pattern 'snr_*.mat' in the
%   current directory (you can change that in the code). Each file should contain
%   a structure named either 'snr_grad' or 'snr_mag'. The loaded structures are
%   compiled into the output variables all_snr_grad and all_snr_mag.
%
% OUTPUT:
%   all_snr_grad - array of structs for grad SNR data
%   all_snr_mag  - array of structs for mag SNR data
%
% -------------------------------------------------------------------------

    % -- (A) Identify all relevant .mat files in your directory -----------
    snrFiles = dir('snr_*.mat'); 
    % If you'd rather load files by hand, you can comment out the line above and 
    % define something like:
    %   snrFiles = {
    %       'snr_grad_tSSSvsRaw_1.mat'
    %       'snr_grad_MEGnet_Vs_tSSS_2.mat'
    %       'snr_mag_MEGnet_Vs_tSSS_2.mat'
    %       'snr_mag_tSSSvsRaw_1.mat'
    %       'snr_mag_tSSSvsRaw_2.mat'
    %   };
    % But the dir() approach will automatically pick up all such files.

    % -- (B) Initialize containers for the SNR data -----------------------
    all_snr_grad = [];  % store SNR grad data in a struct array
    all_snr_mag  = [];  % store SNR mag data in a struct array
    
    % -- (C) Loop through each file and load data -------------------------
    for i = 1:length(snrFiles)
        fprintf('Loading file: %s\n', snrFiles(i).name);
        
        % Load the MAT file into a temporary struct
        tmp = load(snrFiles(i).name);
        
        % Check for presence of snr_grad or snr_mag structure
        if isfield(tmp, 'snr_grad')
            snr_grad_data = tmp.snr_grad;
            
            % Store the data from this file
            all_snr_grad(end+1).filename = snrFiles(i).name;  %#ok<AGROW>
            all_snr_grad(end).sub       = snr_grad_data.sub;
            all_snr_grad(end).sub_run   = snr_grad_data.sub_run;
            all_snr_grad(end).F_all     = snr_grad_data.F_all;
            all_snr_grad(end).snr_value = snr_grad_data.snr_value;

        elseif isfield(tmp, 'snr_mag')
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
    
    % -- (D) Print summary ------------------------------------------------
    fprintf('\nFinished loading:\n');
    fprintf('  SNR grad entries: %d\n', length(all_snr_grad));
    fprintf('  SNR mag entries:  %d\n', length(all_snr_mag));
    
    % -- (E) (Optional) Perform your LI analysis here ---------------------
    % For example:
    %   results_grad = compute_LI_all(all_snr_grad);
    %   results_mag  = compute_LI_all(all_snr_mag);
    % Return them, save them, or display them as needed.
    %
    % If you do your analysis here, you could also return the results:
    %   [all_snr_grad, all_snr_mag, results_grad, results_mag] = ...
    %       ecpfunc_read_noiseSNR();
    %
    % For demonstration well just return the loaded data.

    % (F) End of function - all_snr_grad and all_snr_mag are returned
end
