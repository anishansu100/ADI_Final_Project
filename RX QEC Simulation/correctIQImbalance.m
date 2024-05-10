function [y_corrected] = correctIQImbalance()
    [I_signal, Q_signal] = DirectDownConversionDemo();
    % Apply low-pass filtering
    fbw = 0.2; % bandwidth of the LPF
    N = length(I_signal); % number of samples
    % sample indices, from transmitter
    n = (0:(N-1))';
    % Step 2: Compute βI and βQ (DC offsets)
    beta_I = mean(I_signal);
    beta_Q = mean(Q_signal);
    % Step 3: Remove the DC offsets
    I_error = I_signal - beta_I;
    Q_error = Q_signal - beta_Q;
    %cross_correlationIQ = (alpha*sin(psi) * I_signal.^2);
   
    % Step 4: Compute α (amplitude error) g_q 
    alpha = sqrt(mean(I_error.^2) / mean(Q_error.^2));
    % Step 5: Compute sin(ψ) (phase error) phi_q
    psi = asin((mean(I_error.*Q_error))/ sqrt(mean(I_error.^2) .* mean(Q_error.^2)));
    % Step 7: Compute the correction matrix parameters
    A = 1 / alpha;
    C = -sin(psi) / (alpha * cos(psi));
    D = 1 / cos(psi);
    % Step 8: Apply the correction
    corrected_signal = zeros(2, N);
    corrected_signal(1, :) = A * (I_error);  % Corrected I, first row, second term goes to 0
    corrected_signal(2, :) = C * (I_error) + D * (Q_error);  % Corrected Q, second row
    I_corr = corrected_signal(1, :);
    Q_corr = corrected_signal(2, :); 
    y_corrected = I_corr + 1i*Q_corr; %returning complex signal 
    
    figure
    plot_ft(I_signal + 1i*Q_signal, 'Original Signal');
    plot_ft(y_corrected, 'Corrected Signal');
end

function plot_ft(x, titleStr)
    N = length(x);
    f = (-N/2):(N-1)/2;
    f = f / N; % Normalize frequencies
    
    X = fftshift(fft(x));
    figure;
    plot(f, 20*log10(abs(X)/max(abs(X))));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title(titleStr);
end
