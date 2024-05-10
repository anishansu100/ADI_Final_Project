function plot_file(filename, name)
    Fs = 30.72e6; % Sampling frequency (e.g., 100 kHz
    data_xn = csvread(filename);
    x_n1 = data_xn(:, 1);
    x_n2 = data_xn(:, 2);
    x_n = x_n1 + 1i * x_n2;
    figure
    PlotPsd(x_n, Fs);
    title(name);
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