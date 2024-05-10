function correctIQ_imbalance_t()
    theta_obs = 0; % Observation phase example value (30 degrees converted to radians)
    g_obs = 1; % Gain for the observation channel (example value

    x_n = signal();
    [y_n] = r_to_t_with_error(x_n); % Assuming this function returns x_n and y_n
    [a1, a2] = solve_for_a1_a2(x_n, y_n);
    complex_number = (a1 - a2)/(a1 + a2);
    gq = abs(complex_number);
    theta_q = angle(complex_number);
    % Now calculate 'b' using the values of gq and theta_q
    b = 1/(gq * cos(theta_q)) - 1j * tan(theta_q);

    % Calculate the corrected x_t
    x_t = 1/2 * (1 + b) * x_n + 1/2 * (1 - b) * conj(x_n);

    % Calculate g1 and g2 
    g1 = (1/2) * (1 + gq * exp(1j*theta_q));
    g2 = (1/2) * (1 - gq * exp(1j*theta_q));

    % Calculate y[n] 
    y_new = (1/2) * g_obs * exp(1j*theta_obs) * (g1 * x_t + g2 * conj(x_t));

    % Now plot the corrected signal
    plot_ft(y_new, 'Corrected Signal Spectrum', 100e3); % Assume Fs is 100 kHz as an example
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
