function idea_reciever_WithQE()
    % number of samples to generate
    N = 1e5;
    
    % sample indices
    n = (0:(N-1))';
    
    % Define the carrier frequency and phase
    w_c = 2*pi*0.3; % normalized angular frequency of the carrier
    phi_c = 0; % phase of the carrier in radians
    
    f_x = 0.05;
    % Generate a complex exponential baseband signal
    x_t = exp(1i*2*pi*f_x*n);
    
    % Add some noise to make the plots look more realistic
    x_t = x_t + 10^(-50/20)*sqrt(2)/2*(randn(N,1) + 1i*randn(N,1));
    
    % Create the RF bandpass signal y(t)
    y_t = real(x_t .* exp(1i*(w_c*n + phi_c)) + conj(x_t) .* exp(-1i*(w_c*n + phi_c)));
    
    % Define quadrature error parameters
    g_q = 0.95; % Amplitude imbalance
    phi_q = pi/10; % Phase imbalance in radians
    
    % Calculate g1 and g2 for the I and Q channels with errors
    g1 = 0.5 * (1 + g_q * exp(-1i * phi_q));
    g2 = 0.5 * (1 - g_q * exp(1i * phi_q));
    
    % Down-convert the signal
    y_i_t = real(x_t) .* cos(w_c*n + phi_c); % In-phase component
    y_q_t = -g_q * imag(x_t) .* sin(w_c*n + phi_c); % Quadrature component with error
    
    % Combine in-phase and quadrature components to form complex baseband signal
    y_bb_t = (y_i_t + 1i*y_q_t) .* exp(1i*phi_c);
    
    % Apply complex gains g1 and g2 to I and Q channels
    y_bb_t = real(y_bb_t) * g1 + imag(y_bb_t) * g2;
    
    % Apply low-pass filtering (simple moving average filter)
    lpf_window_size = 100; % Window size of the moving average filter
    y_bb_lpf = filter(ones(lpf_window_size,1)/lpf_window_size, 1, y_bb_t);
    
    % FFT of the filtered signal
    Y_bb_lpf = fftshift(fft(y_bb_lpf));
    
    % Compute frequency vector for plotting
    f = (-N/2):(N-1)/2;
    f = f / N; % Normalize frequencies
    
    % Plot the original and down-converted signal spectra
    figure;
    subplot(3,1,1);
    plot(f, 20*log10(abs(fftshift(fft(x_t)))/N));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title('|X(f)|^2');
    
    subplot(3,1,2);
    plot(f, 20*log10(abs(fftshift(fft(y_bb_t)))/N));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title('|Y_{bb}(f)|^2 Before LPF');
    
    subplot(3,1,3);
    plot(f, 20*log10(abs(Y_bb_lpf)/N));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title('|Y_{bb}(f)|^2 After LPF with Quadrature Error');
end
