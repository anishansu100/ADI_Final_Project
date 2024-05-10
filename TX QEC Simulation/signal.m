function [x_n] = signal()
    % Define parameters
    Fs = 100e3; % Sampling frequency (e.g., 100 kHz)
    N = 1e4; % Number of samples
    t = (0:N-1)'/Fs; % Convert sample indices to time in seconds
    % Baseband I and Q signals 
    x_i = cos(2*pi*1e3*t); % 1 kHz tone for I
    x_q = sin(2*pi*1e3*t); % 1 kHz tone for Q
    % Calculate the complex signal x[n]
    x_n = x_i + 1j*x_q;
end