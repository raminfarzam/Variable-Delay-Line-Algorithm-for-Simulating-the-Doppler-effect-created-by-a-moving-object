function [out_signal] = Variable_Delay(signal , delay, interpolation)

    % Function to introduce variable delay line 
    % (Doppler effect or frequency shift due to sound source movement) 
    % in a signal 'signal' based on possible delay 'delay'
    
    % Important note :
    % if the signal should be delayed with respect to the previous sample
    % the function considers a new sample in between and fill it with liner
    % inter polation
    
    
    % Inputs:
    % - s: input signal (array)
    % - d_sample: delay samples (array of integers indicating delay at each point)
    % - interpolation: default value is linear interpolation can be 
    %                    1- 'linear'
    %                    2- 'spline'
    %                    3- 'pchip' for cubic interpolation
    %                    4- 'next' to add next value of signal
    %                    5- 'previous'  to add the previous value
    % Output:
    % - out_signal: the delayed/rotated version of the input signal
   
    if nargin < 3
        interpolation_type = 'linear';
    end
        
    % Initialize variables
    rpt = 1;                %rpt' keeps track of shifted positions due to delay
    Ns = length(signal);    % Length of input signal
    
    % Preallocate buffers (buff) and indices for NaN handling
    clear buff nanIndices nonNaNIndices
    
    % Loop through each element of the signal
    for i = 1:length(1:Ns)
        buff(i) = signal(i); % Store the input signal in a buffer for manipulation
        
        % Case 1: Handle the first sample
        if i == 1
            % If there's no delay at the first sample, copy it directly
            if delay(i) == 0
                out_signal(i) = buff(i);
                rpt = rpt+1;            % Increment rpt as no delay occurred
                continue;
            else
                % If there is a delay, assign NaN for this position
                out_signal(i) = NaN;
                continue;
            end
        end
        
        % Case 2: No delay between consecutive samples
        if delay(i) == 0 && delay(i-1)== 0
            out_signal(i) = buff(i);        % Copy signal directly
            rpt=rpt+1;
            continue;
        end
        
        % Case 3: Delay has occurred, but check for variations
        if (i - rpt) ~= delay(i)
            if delay(i) > delay(i-1) || delay(i) == delay(i-1)
                % If current delay is same or greater than the previous
                % one(positive delay line)
                out_signal(i)= NaN;         % Assign NaN because of delay
                continue;
            elseif (delay(i) < delay(i-1)) 
                % Increment rpt if the delay decreases (negative delay line)
                rpt = rpt + 1;
            end
        end
         % Case 4: Handle cases where delay matches previous ones
        if (i - rpt) == delay(i) && i>1
            if delay(i)== 0 && delay(i-1) == 1
                % Special case: if current delay is zero but previous was one
                rpt = rpt-1;                    % Adjust rpt accordingly
                out_signal(i) = buff(rpt);      % Get the delayed sample
                rpt = rpt+1;
            else
                rpt = rpt +1;
                out_signal(i) = buff(rpt);      % Assign delayed sample
            end
        end
    end
    % Fill missing NaN values in the signal using linear interpolation
    out_signal = fillmissing(out_signal , interpolation_type);
    % Handle remaining NaNs (if any) using interpolation
    nanIndices = isnan(out_signal);             % Find indices of NaN values
    nonNaNIndices = find(~isnan(out_signal));   % Find indices of valid values
    out_signal(nanIndices) = interp1(nonNaNIndices, out_signal(nonNaNIndices),...
                                     find(nanIndices), interpolation_type);
end

