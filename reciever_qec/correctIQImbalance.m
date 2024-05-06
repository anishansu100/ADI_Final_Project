function [corrected_signal] = correctIQImbalance()
    [I_signal, Q_signal] = DirectDownConversionDemo();
    N = 1e5; % number of samples
    % sample indices, from transmitter
    n = (0:(N-1))';

    % Step 2: Compute βI and βQ (DC offsets)
    beta_I = mean(I_signal);
    beta_Q = mean(Q_signal);

    % Step 3: Remove the DC offsets
    I_error = I_signal - beta_I;
    Q_error = Q_signal - beta_Q;

    %I_error = alpha * I_signal;
    %Q_error = (sin(psi) * I_signal) + (cos(psi) * Q_signal);

    %autocorrelation_I = alpha^2 * I_signal.^2;
    %autocorrelation_Q = I_signal.^2;
    %cross_correlationIQ = (alpha*sin(psi) * I_signal.^2);
   
    % Step 4: Compute α (amplitude error) g_q 
    %alpha = sqrt(2 * I.^2);
    alpha = sqrt(mean(I_error.^2) / mean(Q_error.^2));

    % Step 5: Compute sin(ψ) (phase error) phi_q
    %sin_psi = (2 / alpha) * I .* Q;

    psi = asin((mean(I_error.*Q_error))/ sqrt(mean(I_error.^2) .* mean(Q_error.^2)));

    % Step 6: Compute cos(ψ)
    %cos_psi = sqrt(1 - sin_psi^2);

    % Step 7: Compute the correction matrix parameters
    A = 1 / alpha;
    C = -sin(psi) / (alpha * cos(psi));
    D = 1 / cos(psi);

    % Step 8: Apply the correction
    corrected_signal = zeros(2, N);
    corrected_signal(1, :) = A * (I_signal - I_error);  % Corrected I, first row, second term goes to 0
    corrected_signal(2, :) = C * (I_signal - I_error) + D * (Q_signal - Q_error);  % Corrected Q, second row
end