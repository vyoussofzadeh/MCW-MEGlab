graph TD
    A[Start: Dataset (MEG and fMRI)] --> B[Task/Conditions]
    B --> C1[MEG: Semantic Decision]
    B --> C2[fMRI: Tone decision ]
    
    A --> D[Subjects]
    D --> E1[Patients with Temporal Lobe Epilepsy]
    D --> E2[Healthy Control]
    
    A --> F[Preprocessing]
    F --> G1[MEGnet ICA-based]
    F --> G2[40 Hz low-pass]
    F --> G3[PCA-based ECOG and ECG]
    F --> G4[Rejecting Bad Trials (z-score and kurtosis)]
    F --> G5[SNR Trial Analyses]
    
    A --> H[Source Analysis]
    H --> I1[Time-resolved/evoked dSPM MNE]
    H --> I2[Beta-modulation]
    H --> I3[DICS]
    H --> I4[Hilbert (envelope), induced]
    
    A --> J[Contrasts]
    J --> K1[Active-vs-Baseline]
    J --> K2[Active-vs-Control]
    J --> K3[Active]
    
    A --> L[Compute LI]
    L --> M1[Magnitude]
    L --> M2[Counting]
    L --> M3[Bootstrapping]
    
    A --> N[Atlas and ROIs]
    N --> O1[Atlas: HCP_MMP1]
    N --> O2[ROIs: Angular, Frontal, Temporal, Lateral]
    
    A --> P[Comparisons and Interpretations]
    P --> Q1[LI Concordance & Confusion Matrix]
    P --> Q2[LI Correlation (MEG LIs in time vs. fMRI LI)]
    P --> Q3[Task Performance Scores and Neuropsych/Epilepsy Measures]

    style A fill:#f9f,stroke:#333,stroke-width:2px
    style Q3 fill:#ccf,stroke:#333,stroke-width:2px
