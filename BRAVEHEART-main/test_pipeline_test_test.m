% Define the folder where the DICOM files are located
%ecg_folder = 'C:\\Users\\arthu\\Downloads\\BRAVEHEART-main-20250211T004415Z-001\\BRAVEHEART-main\\AI - ECG Risk Tool';
ecg_folder = "C:\Users\arthu\Downloads\BRAVEHEART-main-20250211T004415Z-001\BRAVEHEART-main\data";
%ecg_folder = "C:\Users\arthu\Downloads\BRAVEHEART-main-20250211T004415Z-001\BRAVEHEART-main\output_new";


% Folder containing the DICOM files
%format = 'DICOM';                 % Set format to 'DICOM' for the DICOM ECG file
format = 'physionet_dat';
output_ext = '.csv';              % Choose '.csv' as the output file extension
output_note = '_patient_';   % Base note for the output file

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

% Get a list of all data files in the folder based on format
disp(['Looking for files in: ' ecg_folder]);
if strcmp(format, 'DICOM')
    input_files = dir(fullfile(ecg_folder, '*.dcm'));
    disp(['Found ' num2str(length(input_files)) ' DICOM files']);
elseif strcmp(format, 'physionet_dat')
    input_files = dir(fullfile(ecg_folder, '*.dat'));
    disp(['Found ' num2str(length(input_files)) ' Physionet DAT files']);
else
    % For other formats, try to get all files
    input_files = dir(fullfile(ecg_folder, '*.*'));
    % Filter out directories and known non-data files
    input_files = input_files(~[input_files.isdir]);
    disp(['Found ' num2str(length(input_files)) ' files']);
end
% Initialize results cell array
all_features = {};

% Call the braveheart_batch function
[results, ecg_raw, filtered_ecg] = braveheart_batch_my_own(ecg_folder, format, output_ext, output_note, parallel_proc, ...
    progressbar, save_figures, save_data, save_annotations, ...
    vcg_calc_flag, lead_morph_flag, vcg_morph_flag);

% Process the result into feature values
result = AnnoResult(results);
[h, a] = result.export_data();

% Number of rows in 'a' corresponds to the number of patients/files processed
num_patients = size(a, 1);
% 
[rms_min_all, rms_var_all] = calculateNoiseMetrics(ecg_raw, filtered_ecg);

% Loop through each patient/file
for patient_idx = 1:num_patients
    % Initialize feature value row
    feature_value = cell(1, 76);
    
    % Get the filename for this patient
    % Since we can't rely on the order of files matching the order of processing,
    % we need to get the filename from the results array if possible
    
    % Try to get the filename from the results data structure if available
    % The exact field depends on your implementation of braveheart_batch_my_own
    % This is a common approach in many batch processing functions
    if isfield(results{patient_idx}, 'filename') || (isstruct(results{patient_idx}) && isfield(results{patient_idx}, 'filename'))
        % If the results structure includes the filename
        [~, filename, ~] = fileparts(results{patient_idx}.filename);
    elseif isfield(ecg_raw{patient_idx}, 'filename') || (isstruct(ecg_raw{patient_idx}) && isfield(ecg_raw{patient_idx}, 'filename'))
        % If the ecg_raw structure includes the filename
        [~, filename, ~] = fileparts(ecg_raw{patient_idx}.filename);
    elseif patient_idx <= length(input_files)
        % Fallback to the input_files array if other methods fail
        % This assumes files are processed in the same order they're listed
        [~, filename, ~] = fileparts(input_files(patient_idx).name);
    else
        % Final fallback if we can't determine the filename
        filename = strcat('patient_', num2str(patient_idx));
    end
    
    % Populate feature values
    feature_value{1} = filename; % Use filename instead of patient_number
    feature_value{2} = "need to label";
    feature_value{3} = "need to obtain";
    feature_value{4} = a{patient_idx, 15};

    ecg_raw_pt = ecg_raw{patient_idx};
    PR_interval = PR_interval_check_PT(ecg_raw_pt.I);

    feature_value{5} = PR_interval;
    feature_value{6} = a{patient_idx, 103+1};
    feature_value{7} = a{patient_idx, 55+1};
    feature_value{8} = a{patient_idx, 33+1};
    feature_value{9} = sqrt(a{patient_idx, 30+1}^2 + a{patient_idx, 31+1}^2);
    feature_value{10} = a{patient_idx, 58+1};
    feature_value{11} = a{patient_idx, 57+1};
    feature_value{12} = a{patient_idx, 324+1};
    feature_value{13} = a{patient_idx, 323+1};
    feature_value{14} = a{patient_idx, 396};
    feature_value{15} = a{patient_idx, 32+1};
    feature_value{16} = a{patient_idx, 32+1} / a{patient_idx, 29+1};
    feature_value{17} = a{patient_idx, 398};
    feature_value{18} = a{patient_idx, 399};
    

    feature_value{19} = str2double(a{patient_idx, 356+1});
    feature_value{20} = str2double(a{patient_idx, 357+1});
    feature_value{21} = str2double(a{patient_idx, 358+1});
    feature_value{22} = a{patient_idx, 24};
    feature_value{23} = a{patient_idx, 350+1};
    feature_value{24} = a{patient_idx, 382+1};
    feature_value{25} = a{patient_idx, 383+1};
    feature_value{26} = a{patient_idx, 371+1};
    feature_value{27} = a{patient_idx, 95+1};
    feature_value{28} = a{patient_idx, 366+1};
    feature_value{29} = a{patient_idx, 368+1};
    feature_value{30} = a{patient_idx, 370+1};
    feature_value{31} = a{patient_idx, 352+1}; 

    temp1 = str2double(a{patient_idx, 337})^2;
    temp2 = str2double(a{patient_idx, 336})^2;

    feature_value{32} = temp1/ temp2;
    feature_value{33} = a{patient_idx, 397};

    temp3 = str2double(a{patient_idx, 339+1})^2;
    temp4 = str2double(a{patient_idx, 338+1})^2;
    feature_value{34} = temp3 / temp4;
    
    temp5 = str2double(a{patient_idx, 394})^2;
    temp6 = str2double(a{patient_idx, 395})^2;
    feature_value{35} = temp6/temp5;

    feature_value{36} = a{patient_idx, 365+1};
    feature_value{37} = a{patient_idx, 367+1};
    feature_value{38} = a{patient_idx, 369+1};
    feature_value{39} = a{patient_idx, 351+1};

    RMS_value = sqrt((1/12)*(a{patient_idx, 180+1}^2 + a{patient_idx, 190+1}^2 + ...
        a{patient_idx, 200+1}^2 + a{patient_idx, 210+1}^2 + a{patient_idx, 220+1}^2 + ...
        a{patient_idx, 230+1}^2 + a{patient_idx, 240+1}^2 + a{patient_idx, 250+1}^2 + ...
        a{patient_idx, 260+1}^2 + a{patient_idx, 270+1}^2 + a{patient_idx, 280+1}^2));
    
    feature_value{40} = RMS_value;
    feature_value{41} = a{patient_idx, 400};
    feature_value{42} = rms_min_all(patient_idx);
    feature_value{43} = rms_var_all(patient_idx);
    feature_value{44} = a{patient_idx, 403};
    feature_value{45} = a{patient_idx, 390+1};
    feature_value{46} = a{patient_idx, 391+1};
    feature_value{47} = a{patient_idx, 393};
    feature_value{48} = a{patient_idx, 196}*1000;
    feature_value{49} = a{patient_idx, 216}*1000;
    feature_value{50} = a{patient_idx, 266}*1000;
    feature_value{51} = a{patient_idx, 385};
    feature_value{52} = a{patient_idx, 386};
    feature_value{53} = a{patient_idx, 387};
    feature_value{54} = a{patient_idx, 388};
    feature_value{55} = str2double(a{patient_idx, 373})*1000;
    feature_value{56} = str2double(a{patient_idx, 374})*1000;
    feature_value{57} = str2double(a{patient_idx, 375})*1000;
    feature_value{58} = str2double(a{patient_idx, 376})*1000;
    feature_value{59} = str2double(a{patient_idx, 377})*1000;
    feature_value{60} = str2double(a{patient_idx, 378})*1000;
    feature_value{61} = str2double(a{patient_idx, 379})*1000;
    feature_value{62} = str2double(a{patient_idx, 380})*1000;
    feature_value{63} = str2double(a{patient_idx, 381})*1000;
    feature_value{64} = str2double(a{patient_idx, 382})*1000;
    feature_value{65} = a{patient_idx, 181}*1000;
    feature_value{66} = a{patient_idx, 191}*1000;
    feature_value{67} = a{patient_idx, 201}*1000;
    feature_value{68} = a{patient_idx, 231}*1000;
    feature_value{69} = a{patient_idx, 221}*1000;
    feature_value{70} = a{patient_idx, 211}*1000;
    feature_value{71} = a{patient_idx, 241}*1000;
    feature_value{72} = a{patient_idx, 251}*1000;
    feature_value{73} = a{patient_idx, 261}*1000;
    feature_value{74} = a{patient_idx, 271}*1000;
    feature_value{75} = a{patient_idx, 281}*1000;
    feature_value{76} = a{patient_idx, 291}*1000;

    % Append to all_features
    all_features = [all_features; feature_value];
end


% Define headers - change the first column name to "Filename" instead of "patient_1"
features_header = {"Filename","OMI", "Age", "HR", "PR", "QRSd", "mfpQRSaxis", ...
    "fpTaxis", "txzQRSaxis", "QRSTangle", "mQRSTangle", "TCRTangle", "TCRT", ...
    "TampInfl1", "Tamp", "TrelAmp", "fpTinfl1Axis", ...
    "fpTinfl1Mag", "PCA1", "PCA2", "PCA3", "spatialTaxisDev", ...
    "TMD", "TMDpre", "TMDpost", "JTc", "TpTe", "pctRNDPV", "pctJNDPV", ...
    "pctSTNDPV", "pctTNDPV", "QRS_PCAratio", "J_PCAratio", "STT_PCAratio", ...
    "T_PCAratio", "RNDPV", "JNDPV", "STNDPV", "TNDPV", "pcaTamp", "Tasym", ...
    "rmsMin", "rmsVar", "MIsz", "latConcaveAmp", "antConcaveAmp", "qdur_aVF", ...
    "ramp_III", "ramp_aVL", "ramp_V4", "vat_II", "vat_III", "vat_V2", ...
    "vat_V4", "st80_I", "st80_III", "st80_aVL", "st80_aVF", "st80_V1", ...
    "st80_V2", "st80_V3", "st80_V4", "st80_V5", "st80_V6", "tamp_I", ...
    "tamp_II", "tamp_III", "tamp_aVR", "tamp_aVL", "tamp_aVF", "tamp_V1", ...
    "tamp_V2", "tamp_V3", "tamp_V4", "tamp_V5", "tamp_V6"};

% Save the results to a CSV file
outname = strcat('processed_ecgs_all', output_ext);  % Output file name
outputfilename = fullfile(ecg_folder, outname);
writecell([features_header; all_features], outputfilename);