function [y_new] = r_to_t_with_error(x_n, title, gq, theta_q)
    % Define parameters
    Fs = 30.72e6; % Sampling frequency (e.g., 100 kHz
    f_tx = 2.005e9; % Transmission frequency 
    f_rx = 2e9; % Reception frequency 

    % Define the quadrature error parameters
    theta_obs = 0; % Observation phase example value (30 degrees converted to radians)
    g_obs = 1; % Gain for the observation channel (example value)

    % Calculate g1 and g2 
    g1 = (1/2) * (1 + gq * exp(1j*theta_q));
    g2 = (1/2) * (1 - gq * exp(1j*theta_q));
    
    % Calculate y[n] 
    y_n = (1/2) * g_obs * exp(1j*theta_obs) * (g1 * x_n + g2 * conj(x_n));
    n = (0:(length(x_n)-1))'; % Time index array
    freq_shift = exp(1j * 2 * pi * ((f_tx - f_rx)/Fs) * n );
    y_n = y_n .* freq_shift;
    y_new = circshift(y_n, randi([0, length(y_n)-1], 1, 1));
    % Plot the spectrum of y[n]
    plot_ft(y_new, title, Fs);
end

function plot_ft(x, titleStr, Fs)
    N = length(x);
    f = (-N/2):(N/2-1);
    f = f * Fs / N; % Convert to actual frequency in Hz
    
    X = fftshift(fft(x));
    figure;
    plot(f, 20*log10(abs(X)/max(abs(X))));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (Hz)');
    title(titleStr);
end
