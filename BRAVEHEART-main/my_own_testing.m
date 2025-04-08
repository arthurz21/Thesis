% Define the folder where the DICOM file is located
ecg_folder = "C:\Users\arthu\Downloads\BRAVEHEART-main-20250113T060747Z-001\BRAVEHEART-main\AI - ECG Risk Tool";

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  EDIT BELOW HERE AT YOUR OWN RISK!!!                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up any errant progress bars
all_waitbars = findall(0,'type','figure','tag','TMWWaitbar');
delete(all_waitbars);

% If deployed as executable, will need to pull the configuration parameters
% from external file batch_settings.csv
% This is DISABLED if run from compiled GUI as pass in extra parameter to
% disable
if isdeployed && ~exist('disable_read','var')  % This only is true if running from compiled command line
    [format, output_ext, output_note, parallel_proc, progressbar, ...
        save_figures, save_data, save_annotations, ...
        vcg_calc_flag, lead_morph_flag, vcg_morph_flag] = read_batch_settings(fullfile(getcurrentdir(),'batch_settings.csv'));
end

% To deal with progress bar issues with parallel computing will set the
% parallel computing flag 'parallel_proc' = 0 if the user does not have the
% parallel computing toolbox.  This isn't technically necessary for calculations 
% due to the fact that parfor will run as a regular for loop if the toolbox is 
% not installed, but disabling the flag makes dealing with the progress bars
% sigificantly easier as its require calling commands that do not exist
% without the toolbox.  This really is the only way to allow both parallel
% and serial waitbars to work

if isempty( ver('parallel')) 
    parallel_proc = 0;
end

tic     % Start timer

% Extension of annotation files
anno_ext = '.anno';     

% Get file extension
ecg_source_hash = read_format_csv(fullfile(getcurrentdir(),'ecg_formats.csv'));
[~, source_ext, ~] = ecg_source_string(format, ecg_source_hash);

% Date and time formatting for saving
formatDate = 'mmddyyyy';
formatTime = 'HHMMSS';
time_stamp = strcat('_',datestr(now,formatDate),'_',datestr(now,formatTime));

% Load folder of ECGs:

% If auto source extension, use source_ext helper function to get the
% expected extension for the file based on source_str
if strcmp(source_ext, 'auto')
   source_ext = get_source_ext(format);   
end

% Progress bar stuff
if progressbar
    H = waitbar(0, 'Please wait ... Initializing','Name','Processing...');
    
    % Cant create DataQueue if don't have parallel computing toolbox
    if parallel_proc
        D = parallel.pool.DataQueue;
        afterEach(D, @UpdateWaitbar);
    else 
        clear D
        D = 0;
    end
else 
    D = 0;
    H = 0;
end

% Folder file stucture
file_list_struct = dir(fullfile(ecg_folder, strcat('*',source_ext)));
name_list_struct = {file_list_struct.name}';
name_list_struct(ismember(name_list_struct, {'.', '..'})) = [];

% Create directory names
orig_directory = file_list_struct.folder;
fig_full_directory = fullfile(orig_directory, strcat('figures', time_stamp));
checkfig_directory = fullfile(orig_directory, strcat('figures', time_stamp,'_check_ecg_list'));
annotation_directory = fullfile(orig_directory, strcat('annotations', time_stamp));
data_directory = fullfile(orig_directory, strcat('data', time_stamp));

% Load ECGs from directory
num_files = length(name_list_struct);
file_list = cell(num_files, 1);
for i = 1:num_files
    file_list{i} = fullfile(orig_directory,name_list_struct(i));
end

% Create directories based on user input
% Create directory if needed to save figures
if save_figures  
	mkdir(fig_full_directory);
	mkdir(checkfig_directory);
end

% Create directory if needed to save annotation files
if save_annotations
	mkdir(annotation_directory); %#ok<*UNRCH>
end

% Create directory if needed to save data
if save_data
	mkdir(data_directory); %#ok<*UNRCH>
end

% Create flags structure
flags = struct;
flags.vcg_calc_flag = vcg_calc_flag;
flags.lead_morph_flag = lead_morph_flag;
flags.vcg_morph_flag = vcg_morph_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize quality and error counts
qindex = zeros(num_files,1);
eindex = zeros(num_files,1);

% If have passed in annoparams will use them here.  
% If have not passed in annoparams, then will use the default settings that
% are pulled from Annoparams.m or (Annoparams.csv if deployed)
if exist('ap','var')
    default_ap = ap;
else
    default_ap = Annoparams();
end

% Load Qualparams (will load Qualparams.csv if deployed)
qp = Qualparams();

% If want to change median beat annotater (not recomended!)
% default_ap.median_reanno_method = 'Std';

% Preallocate output
results = cell(num_files, 1);
quality_mat = cell(num_files, 1);
quality_full = cell(num_files,1);
error_ecgs = cell(num_files,1);
error_msgs = cell(num_files,1);
exportsigs = strings(num_files,1);

% Open parallel pool if specified use of parallel computing
if parallel_proc
    workers = Inf;
else
    workers = 0;
end

P = 1;      % Parallel iterating counter
ecg = ECG12(char(file_list{1}), format);

num_files = 1;
i=1;


note = output_note;

% Parse out filename and extension -- don't want extension in file/figure names
[~, basename, fext] = fileparts(file_list{i});
basename = char(basename);
fext = char(fext);

try   % To deal with if ECG wont load due to some error in formatting
	
	% Load ECG from file and source_str
    ecg = ECG12(char(file_list{i}), format)

	% Load custom annotation / anno parameters if any
    % annotation/anno params file has to have the same basename as the
    % ECG, and end in whatever anno_ext is set above, and be located
    % in the same directory as the ecg file
    % ex: ecg123.xml & ecg123.csv
	apfile = strcat(fullfile(orig_directory, basename), anno_ext);
	
    % If the apfile exists, parse out the beats and anno parameters
	if exist(char(apfile), 'file')
		[ap, beats] = Annoparams(apfile);
		ap.debug = false;
        ap.outlier_removal = 0;      % Dont want to remove PVCs and outliers automatically, because this has been done manually already
        ap.pvc_removal = 0;

		if isa(beats, 'Beats')
			note = 'manual';        % Lets you know this ECG had manual editing performed
            
            % If apfile includes beats, then pass beats & anno parameters into batch_calc
			
            [hr, num_initial_beats, ~, quality, corr, medianvcg1, beatsig_vcg, median_12L, beatsig_12L, medianbeat, beat_stats, ecg_raw, vcg_raw, filtered_ecg, filtered_vcg, noise, sumfig] = ...
            batch_calc(ecg, beats, [], [], [], [], ap, qp, save_figures, basename, []);
            %median_12L

            [geh, lead_morph, vcg_morph] = module_output(median_12L, medianvcg1, medianbeat, ap, flags);

        else
            
            % If apfile only includes anno parameters and no beats, pass the annoparameters into batch_calc
            
            [hr, num_initial_beats, beats, quality, corr, medianvcg1, beatsig_vcg, median_12L, beatsig_12L, medianbeat, beat_stats, ecg_raw, vcg_raw, filtered_ecg, filtered_vcg, noise, sumfig] = ...
            batch_calc(ecg, [], [], [], [], [], ap, qp, save_figures, basename, []);
            %median_12L
			    
        [geh, lead_morph, vcg_morph] = module_output(median_12L, medianvcg1, medianbeat, ap, flags);
        
        end
        
    else        % Are not suopplying an apfile with anno parameters and/or beats
		ap = default_ap;
		beats = [];
        
        % Pass ecg and anno parameters into batch_calc
        [hr, num_initial_beats, beats, quality, corr, medianvcg1, beatsig_vcg, median_12L, beatsig_12L, medianbeat, beat_stats, ecg_raw, vcg_raw, filtered_ecg, filtered_vcg, noise, sumfig] = ...
            batch_calc(ecg, beats, [], [], [], [], ap, qp, save_figures, basename, []);

        %median_beats_test = median_12L
        
        [geh, lead_morph, vcg_morph] = module_output(median_12L, medianvcg1, medianbeat, ap, flags);

    end
    

    QRST_pts = [medianbeat.Q medianbeat.QRS medianbeat.S medianbeat.T medianbeat.Tend];
    R_wave_limits = index_finder_vcg(medianvcg1, QRST_pts);
    [RWR_abs, RWR_rel] = svd_rwave(median_12L, R_wave_limits);
    %now we work on ST_J wave Residuum. 

    %[R_wave_limits, RWR_abs, RWR_rel] = svd_rwave_st_wave(median_12L, QRST_pts);
    %[vat_II, vat_III, vat_V2, vat_V4] = find_VAT(median_12L, QRST_pts);

    %we now check for 
    
    
	% At this point, all data is generated - need to organize for export
    vcg_morph;
    % Generate row of results for output file
	results{i} = AnnoResult(strcat(basename,fext), note, format, ap, ecg, hr, num_initial_beats, beats, beat_stats, corr, noise, quality.prob_value, quality.missing_lead, geh, lead_morph, vcg_morph);
	
    % Signal quality object
	quality_mat{i} = quality;
	
    % Figure file names/path for saving
	fig_filename = char(strcat(basename,'.png'));
	fig_full_path = char(fullfile(fig_full_directory, fig_filename));
	
    % If some quality object was flagged, will add the ECG/figure to
    % the check ECG list
	if quality.counter() >= 1
		quality_full{i} = [strcat(basename,fext) num2cell(quality_mat{i}.vector())];
		qindex(i) = 1;  % Increment counter for number of ECGs that were flagged for quality issues
		
		% Save figure in separate directory for easy viewing
		if save_figures; saveas(sumfig,fullfile(checkfig_directory, fig_filename)); end	
	end
	
	% Save figure
	if save_figures
		saveas(sumfig,fig_full_path);
		close(sumfig)
	end
	
	if save_data

    	sig = struct( ...
            'filename',file_list{i}, ...
            'annoparams', ap, ...
            'ecg_raw', ecg_raw, ...
            'ecg_filtered', filtered_ecg, ...
            'vcg_raw', vcg_raw, ...
            'vcg_filtered', filtered_vcg, ...
            'beats', beats, ...
            'beat_stats', beat_stats, ...
            'geh', geh, ...
            'median_vcg', medianvcg1, ...
            'median_12L', median_12L, ...
            'medianbeat', medianbeat, ...
            'beats_median_vcg', beatsig_vcg, ...
            'beats_median_12L', beatsig_12L, ...
            'lead_morph', lead_morph, ...
            'vcg_morph', vcg_morph, ...
            'quality', quality);
		
		save_struct_parfor(data_directory, basename, sig)
		
	end
	
	if save_annotations
		annotation_filename = char(strcat(basename, anno_ext));
		annotation_full_path = char(fullfile(annotation_directory, annotation_filename));
		ap.to_file(beats, annotation_full_path);
	end

	if isnan(geh.svg_x) && vcg_calc_flag == 1
	  error_ecgs{i} = basename; 
	  error_msgs{i} = 'VCG Calc calculation failed';
	  eindex(i) = 1;	
	end

	if isnan(lead_morph.L1_r_wave) && lead_morph_flag == 1
	  error_ecgs{i} = strcat(basename,fext); 
	  error_msgs{i} = 'Lead Morphology calculation failed';
	  eindex(i) = 1;	
    end
    
    if isnan(vcg_morph.TCRT) && vcg_morph_flag == 1
	  error_ecgs{i} = strcat(basename,fext); 
	  error_msgs{i} = 'VCG Morphology calculation failed';
	  eindex(i) = 1;	
	end
    
catch ME  % error catching as part of try loop
	error_ecgs{i} = strcat(basename,fext); 
	error_msgs{i} = ME.message;
	eindex(i) = 1;	
end	

% Update progress bars
if progressbar
    if parallel_proc
        send(D, i);          % Send iteration counter to DataQueue
    else
        k = num_files - i+1  % Turns out if there are 0 workers parfor seems to go in reverse order!
        waitbar(k/num_files, H, sprintf('Processed %i out of %i Total ECGs (%i%%)',k,num_files,round(100*(k/num_files))),'Name','Processing...');
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write all data to Excel or CSV 
outname = strcat('processed_ecgs',time_stamp,output_ext);  % Excel spreadhseet or CSV filename to save results
outputfilename = fullfile(orig_directory,outname);

% Export to .xlsx or .csv
result = AnnoResult(results);

[h, a] = result.export_data();
writecell( [h ; a], outputfilename);



% Export to .xlsx or .csv

% Create file for poor quality ecgs to be verfied
if sum(qindex) > 0
    
	% Write list of ecgs that need quality verification (check_ecg_list.xlsx)
	check_ecg_filename = fullfile(orig_directory,'check_ecg_list');
	check_ecg_filename = strcat(check_ecg_filename, time_stamp,output_ext);
	check_ecg_list_header{1}=[{'filename'}, {'qt_error'},{'qrs_error'},{'tpeakqt_error'},...
        {'tmag_error'},{'HR_error'},{'num_beats_error'},{'num_removed_beats_error'},...
        {'median_sig_corr_error'},{'baseline'},{'missing_lead'},{'hf_noise'},{'lf_noise'},...
        {'probability'}, {'NNet_prob_flag'}, {'NNet_NaN'}];
	quality_full = quality_full(~cellfun('isempty',quality_full));
	check_ecg_excel_data = [check_ecg_list_header ; quality_full];
	check_ecg_excel_data = vertcat(check_ecg_excel_data{:});
	writecell( check_ecg_excel_data, check_ecg_filename);

else        % If no quality flags can delete the check_ecg_list directory
    if save_figures
        rmdir(checkfig_directory);
    end
end

if sum(eindex) > 0
    
	% Create list of ecgs that could not load/process due to errors
	error_filename = fullfile(orig_directory,'errors');
	error_filename = strcat(error_filename, time_stamp,output_ext);
	error_ecgs = error_ecgs(~cellfun('isempty',error_ecgs));
	error_msgs = error_msgs(~cellfun('isempty',error_msgs));
	error = [error_ecgs error_msgs];
	writematrix(string(error), error_filename);
end


time = toc;     % End timer

% Calculate duration of batch run
[~,~,~,H,M,S] = datevec(time/(24*60*60));
time_str = '';
if H > 0
   time_str = strcat(num2str(H),{' '},'Hours',{' '});
end
if M > 0
   time_str = strcat(time_str,{''},num2str(M),{' '},'Minutes',{' '});
end
if S > 0
   time_str = strcat(time_str,{''},num2str(round(S,2)),{' '},'Seconds');
end

% Clean up any errant progress bars
all_waitbars = findall(0,'type','figure','tag','TMWWaitbar');
delete(all_waitbars);

if progressbar
    mb = msgbox({sprintf('%i Total ECGs Processed',num_files);sprintf('%i Need Quality Verification',sum(qindex));sprintf('%i Errors',sum(eindex));sprintf('Total time = %s',string(time_str))},'Complete');
end


% Command line text to indicate processing is complete
fprintf('\n\nECG PROCESSING COMPLETE! \n\n%i ECGs Processed \n%i ECGs Need Quality Verification \n%i Errors\nTotal time = %s\n\n',i, sum(qindex), sum(eindex), string(time_str));

% Function to update the parallel waitbar
function UpdateWaitbar(~)
    waitbar(P/num_files, H, sprintf('Processed %i out of %i Total ECGs (%i%%)',P,num_files,round(100*(P/num_files))));
    P = P + 1;
end


