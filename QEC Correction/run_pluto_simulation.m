function run_pluto_simulation()
    Fs = 30.72e6; % Sampling frequency (e.g., 100 kHz
    f_tx = 2.005e9; % Transmission frequency 
    f_rx = 2e9; % Reception frequency 
    g_obs = 1;
    theta_obs = 0;
    % set to 1 for simulation
    if(1)
       [y_n, gq, theta_q] = correctIQ_imbalance_t("dummy.csv", 0, 1);
    else
        [y_n, gq, theta_q] = correctIQ_imbalance_t("dummy.csv", 1, 0);
        [corrected_signal] = correctIQImbalance(real(y_n), imag(y_n));
    end
   %% Trial of  RX qec
   % [corrected_signal] = correctIQImbalance(real(y_n), imag(y_n));
   % y_i = corrected_signal(:, 1);
   % y_q = corrected_signal(:, 2);
   % result = y_i + 1i * y_q;
   % figure
   % PlotPsd(result, Fs);
end

function PlotPsd(x, fs)

    if (nargin < 2)
        fs = 1;
    end
    
    [pxx, f] = GetPsd(x, fs);
    
    plot(f/1e6, 10*log10(pxx));
    xlabel('Frequency (MHz)');
    ylabel('PSD (dBFS/Hz)');
    
end

function [pxx, f] = GetPsd(x, fs)

    if (nargin < 2)
        fs = 1;
    end

    N = 2^floor(log2(length(x)/16));
    M = N/4;
    L = N*4;
    w = blackman(N);

    [pxx, f] = pwelch(x, w, M, L, fs, 'centered');
    
end