function correctIQ_imbalance_t()
    Fs = 30.72e6; % Sampling frequency (e.g., 100 kHz
    f_tx = 2.005e9; % Transmission frequency 
    f_rx = 2e9; % Reception frequency 

    x_n = GenerateCalWaveform();
    [y_n] = r_to_t_with_error(x_n); % Assuming this function returns x_n and y_n
    y_aligned = Alignment(x_n, y_n, f_tx, f_rx, Fs);

    [a1, a2] = solve_for_a1_a2(x_n, y_aligned);
    complex_number = (a1 - a2)/(a1 + a2);
    gq = abs(complex_number);
    theta_q = angle(complex_number);
    % Now calculate 'b' using the values of gq and theta_q
    b = 1/(gq * cos(theta_q)) - 1j * tan(theta_q);

    % Calculate the corrected x_t
    x_t = 1/2 * (1 + b) * x_n + 1/2 * (1 - b) * conj(x_n);

    [y_new] = r_to_t_with_error(x_t);
end



function [a1, a2] = solve_for_a1_a2(x_n, y_n)
    % Compute the summation terms
    S_xx_star = mean(x_n .* conj(x_n));
    S_yx_star = mean(y_n .* conj(x_n));
    S_yx = mean(y_n .* x_n);
    S_xx = mean(x_n .* x_n);

    % Form the matrix and vector for solving the linear equations
    A = [S_xx_star, conj(S_xx); 
         S_xx, S_xx_star];
    b = [S_yx_star; S_yx];

    % Solve for a1 and a2
    a = A \ b; 
    
    % Extract a1 and a2 from the solution vector
    a1 = a(1);
    a2 = a(2);
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
