% Define the folder where the DICOM file is located
ecg_folder = '/home/arthurz21/Downloads/Thesis/BRAVEHEART-main/AI - ECG Risk Tool';
%ecg_folder = "C:\Users\arthu\Downloads\OneDrive_2024-11-18\AI - ECG Risk Tool_full";
% Folder containing the DICOM file(s)

format = 'DICOM';                 % Set format to 'DICOM' for the DICOM ECG file
output_ext = '.csv';              % Choose '.csv' as the output file extension
output_note = '_patient_1';        % Optional note added to the output file

% Processing flags
parallel_proc = 1;                % Use parallel processing (1 for yes, 0 for no)
progressbar = 1;                  % Display a progress bar (1 for yes, 0 for no)
save_figures = 0;                 % Do not save figures
save_data = 0;                    % Save data file with calculated metrics
save_annotations = 0;             % Do not save annotations

% Analysis flags to enable TCRT calculation
vcg_calc_flag = 1;                % Process VCG calculations
lead_morph_flag = 1;              % Process lead morphology (if needed)
vcg_morph_flag = 1;               % Process VCG morphology, necessary for TCRT

% Call the braveheart_batch function
results = braveheart_batch_my_own(ecg_folder, format, output_ext, output_note, parallel_proc, ...
    progressbar, save_figures, save_data, save_annotations, ...
    vcg_calc_flag, lead_morph_flag, vcg_morph_flag);

outname = strcat('processed_ecgs',output_note,output_ext);  % Excel spreadhseet or CSV filename to save results
outputfilename = fullfile(ecg_folder,outname);

% Export to .xlsx or .csv
result = AnnoResult(results);
[h, a] = result.export_data();
features_header = {"patient_1","OMI", "Age", "HR", "PR", "QRSd", "mfpQRSaxis"...
    "fpTaxis", "txzQRSaxis", "QRSTangle", "mQRSTangle", "TCRTangle", "TCRT"...
    "TampInfl1", "Tamp", "TrelAmp","fpTinfl1Axis"...
    "fpTinfl1Mag", "PCA1", "PCA2", "PCA3", "spatialTaxisDev"...
    "TMD", "TMDpre","TMDpost","JTc","TpTe","pctRNDPV", "pctJNDPV"...
    "pctSTNDPV", "pctTNDPV", "QRS_PCAratio", "J_PCAratio", "STT_PCAratio"...
    "T_PCAratio", "RNDPV", "JNDPV", "STNDPV", "TNDPV","pcaTamp", "Tasym"...
    "rmsMin", "rmsVar", "MIsz", "latConcaveAmp", "antConcaveAmp", "qdur_aVF"...
    "ramp_III", "ramp_aVL", "ramp_V4", "vat_II", "vat_III", "vat_V2"...
    "vat_V4", "st80_I", "st80_III", "st80_aVL", "st80_aVF", "st80_V1"...
    "st80_V2", "st80_V3","st80_V4", "st80_V5", "st80_V6", "tamp_I"...
    "tamp_II", "tamp_III", "tamp_aVR", "tamp_aVL", "tamp_aVF", "tamp_V1"...
    "tamp_V2", "tamp_V3", "tamp_V4", "tamp_V5", "tamp_V6"};

feature_value = cell(1,76);
feature_value{1} = "patient_1";
feature_value{2} = "need to label";
feature_value{3} = "need to obtain";
feature_value{4} = a{15};
feature_value{5} = "need to compute, use paper";
feature_value{6} = a{103+1};
feature_value{7} = a{55+1};
feature_value{8} = a{33+1};
feature_value{9} = sqrt(a{30+1}^2 + a{31+1}^2);
feature_value{10} = a{58+1};
feature_value{11} = a{57+1};
feature_value{12} = a{324+1};
feature_value{13} = a{323+1};
feature_value{14} = "";
feature_value{15} = a{32+1};
feature_value{16} = a{32+1}/a{29+1};
feature_value{17} = "";
feature_value{18} = "";
feature_value{19} = str2double(a{356+1})^2;
feature_value{20} = str2double(a{357+1})^2;
feature_value{21} = str2double(a{358+1})^2;
feature_value{22} = a{24};
feature_value{23} = a{350+1};
feature_value{24} = a{382+1};
feature_value{25} = a{383+1};
feature_value{26} = a{371+1};
feature_value{27} = a{95+1};
feature_value{28} = a{366+1};
feature_value{29} = a{368+1};
feature_value{30} = a{370+1};
feature_value{31} = a{352+1}; 
feature_value{32} = (a{336+1}.^2)/(a{335+1}.^2);
feature_value{33} = "";
feature_value{34} = (a{339+1}.^2)/(a{338+1}.^2);
feature_value{35} = "";
feature_value{36} = a{365+1};
feature_value{37} = a{367+1};
feature_value{38} = a{369+1};
feature_value{39} = a{351+1};

RMS_value = sqrt((1/12)*(a{180+1}^2+a{190+1}^2+a{200+1}^2+a{210+1}^2+a{220+1}^2+a{230+1}^2 ...
+a{240+1}^2+a{250+1}^2+a{260+1}^2+a{270+1}^2+a{280+1}^2));

feature_value{40} = RMS_value;
feature_value{41} = "";
feature_value{42} = a{17};
feature_value{43} = "";
feature_value{44} = "";
feature_value{45} = a{390+1};
feature_value{46} = a{391+1};
feature_value{47} = a{393};
feature_value{48} = a{195+1};
feature_value{49} = a{215+1};
feature_value{50} = a{265+1};
feature_value{51} = a{384+1};
feature_value{52} = a{385+1};
feature_value{53} = a{386+1};
feature_value{54} = a{387+1};
feature_value{55} = a{372+1};
feature_value{56} = a{373+1};
feature_value{57} = a{374+1};
feature_value{58} = a{375+1};
feature_value{59} = a{376+1};
feature_value{60} = a{377+1};
feature_value{61} = a{378+1};
feature_value{62} = a{379+1};
feature_value{63} = a{380+1};
feature_value{64} = a{381+1};
feature_value{65} = a{180+1};
feature_value{66} = a{190+1};
feature_value{67} = a{200+1};
feature_value{68} = a{230+1};
feature_value{69} = a{220+1};
feature_value{70} = a{210+1};
feature_value{71} = a{240+1};
feature_value{72} = a{250+1};
feature_value{73} = a{260+1};
feature_value{74} = a{270+1};
feature_value{75} = a{280+1};
feature_value{76} = a{290+1};



writecell( [features_header ; feature_value], outputfilename);