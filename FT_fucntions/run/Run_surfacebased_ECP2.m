switch method
    case 1
        %% SPM surface-based
%         mridir = fullfile(indir,subj,'Anatomy/mri');
        mridir = outputmridir;
        mripfile = fullfile(mridir,'T1.nii');        
        Run_spm
    case 2
        %% BrainStorm export preprocessed ft-IC
        Run_bs_SM2
    case 3
        %% Other implementations
        vy_surfacebasedsource3_ECP
        vy_surfacebasedsource_dics
end