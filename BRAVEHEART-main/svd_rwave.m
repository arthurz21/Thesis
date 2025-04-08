function [RWR_abs, RWR_rel] = svd_rwave(ecg, R_wave_limits, fidpts)

%start_pt = R_wave_limits.avg_start;
start_pt = fidpts(1);
%end_pt = R_wave_limits.avg_end;
end_pt = fidpts(3);
if end_pt < start_pt +8
    end_pt = start_pt+8;
end 
% Cut out R wave

X = [ecg.I(start_pt:end_pt) ecg.II(start_pt:end_pt) ...
    ecg.V1(start_pt:end_pt) ecg.V2(start_pt:end_pt) ...
    ecg.V3(start_pt:end_pt) ecg.V4(start_pt:end_pt) ...
    ecg.V5(start_pt:end_pt) ecg.V6(start_pt:end_pt)]';

% DONT Subtract out the centroid to avoid offsets since already zerod data

% SVD - Only care about left singular vectors U and singular values, but
% easier to extract singular values using second line of code than to get
% the entire matrix S with all the zeros etc.
s=svd(X);
eigmat = s.^2;

% TWR is sum of 4th through nth (8th in this case) eigenalues
RWR_abs = sum(eigmat(4:8));

% Relative TWR is the percentage of the whole
RWR_rel = 100*RWR_abs / sum(eigmat);
end