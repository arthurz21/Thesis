function t_axis_deviation = computeTAxisDeviation(XT_area, YT_area, ZT_area)
    % Inputs:
    % - XT_area: Mean vector area for the T wave in X direction
    % - YT_area: Mean vector area for the T wave in Y direction
    % - ZT_area: Mean vector area for the T wave in Z direction
    % Output:
    % - t_axis_deviation: Angular deviation of the T-axis in degrees

    % Construct the spatial T vector components
    svg_t_x = XT_area;
    svg_t_y = YT_area;
    svg_t_z = ZT_area;

    % Combine into the T vector
    T_vector = [svg_t_x; svg_t_y; svg_t_z];

    % Define the reference vector (e.g., normal T-axis direction)
    reference_vector = [1/sqrt(3); 1/sqrt(3); -1/sqrt(3)];  % Approximation

    % Normalize the T vector and reference vector
    T_unit = T_vector / norm(T_vector);
    ref_unit = reference_vector / norm(reference_vector);

    % Compute the T-axis deviation
    t_axis_deviation = acosd(dot(T_unit, ref_unit));

    % Display results
    %fprintf('T-axis deviation: %.2f degrees\n', t_axis_deviation)
end
