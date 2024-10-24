function mse_result = calculate_mse(signal1, signal2, normalize)
    % Check if the normalization flag is provided
    if nargin < 3
        normalize = false; % Default to false if not provided
    end

    % Normalize signals if specified
    if normalize
        normalized_signal1 = signal1 / norm(signal1);
        normalized_signal2 = signal2 / norm(signal2);
    else
        normalized_signal1 = signal1;
        normalized_signal2 = signal2;
    end

    % Find the common length
    common_length = min(length(normalized_signal1), length(normalized_signal2));

    % Trim or resize normalized signals to the common length
    trimmed_normalized_signal1 = normalized_signal1(1:common_length);
    trimmed_normalized_signal2 = normalized_signal2(1:common_length);

    % Calculate MSE
    mse_result = sum((trimmed_normalized_signal1 - trimmed_normalized_signal2).^2) / common_length;
end
