data = load('JS00001.mat');
% If you specifically need to maintain 12×5000 format in the output
% (note: this is not standard for WFDB)
signals = data.val; % Your 12×5000 data
fs = 500; % Your sampling frequency
sigNames = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
units = cell(1, 12);
[units{:}] = deal('mV');

% Write directly to binary file
fid = fopen('11.dat', 'wb');
for i = 1:size(signals, 1) % For each channel/lead
    fwrite(fid, signals(i,:), 'int16');
end
fclose(fid);