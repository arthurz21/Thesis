%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRAVEHEART - Open source software for electrocardiographic and vectorcardiographic analysis
% svd_8lead.m -- SVD of 8 leads of an ECG
% Copyright 2016-2024 Hans F. Stabenau and Jonathan W. Waks
% 
% Source code/executables: https://github.com/BIVectors/BRAVEHEART
% Contact: braveheart.ecg@gmail.com
% 
% BRAVEHEART is free software: you can redistribute it and/or modify it under the terms of the GNU 
% General Public License as published by the Free Software Foundation, either version 3 of the License, 
% or (at your option) any later version.
%
% BRAVEHEART is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. 
% If not, see <https://www.gnu.org/licenses/>.
%
% This software is for research purposes only and is not intended to diagnose or treat any disease.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s1, s2, s3] = my_own_svd_8lead(ecg)
    X = [ecg.I ecg.II ecg.V1 ecg.V2 ecg.V3 ecg.V4 ecg.V5 ecg.V6];
    
    % Center the data (subtract mean)
    X_centered = X - mean(X);
    
    % Compute covariance matrix
    C = (X_centered' * X_centered) / (size(X_centered, 1) - 1);
    
    % Compute eigenvectors and eigenvalues of covariance matrix
    [V, D] = eig(C);
    
    % Sort eigenvalues (descending)
    [eigenvals, idx] = sort(diag(D), 'descend');
    
    % Calculate percentage of variance explained by each component
    total_variance = sum(eigenvals);
    percent_var = eigenvals / total_variance * 100;
    
    % Return the percentage of variance explained by first 3 components
    s1 = percent_var(1);
    s2 = percent_var(2);
    s3 = percent_var(3);
end