function rx_simulation()
    [corrected_signal] = correctIQImbalance();
    final_signal = corrected_signal(:, 1) + 1j * corrected_signal(:, 2);
    plot_ft(final_signal, "");
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