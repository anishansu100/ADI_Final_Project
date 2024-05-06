% x    [N,1] Tx input signal
% y    [N,1] Rx output signal
% ftx  [1,1] Tx LO frequency
% frx  [1,1] Rx LO frequency
% fs   [1,1] Sampling rate for 'x' and 'y'
% z    [1,1] Output after aligning 'y' to 'x'
function z = Alignment(x, y, ftx, frx, fs)

    N = length(y);

    % frequency-alignment
    z1 = y(:) .* exp(-1i*2*pi*(ftx-frx)/fs*(1:N)');
    
    % coarse time-alignment
    D = 20;  % works well for 500kHz-wide signal sampled at 30.72MHz
    m = (0:(N-1))';
    zz = zeros(N,1);
    for i = 1:D:length(m)
        z2 = circshift(z1, m(i));
        zz(i) = z2' * (x(:));
    end    
    [~,i] = max(abs(zz));
    z2 = circshift(z1, m(i));
    
    % fine time-alignment
    L = 100;
    [zzz, mm] = xcorr(x, z2, L);    
    [~,j] = max(abs(zzz));
    z = circshift(z2, mm(j));
    
    if (0)
        figure;
        set(gcf, 'WindowStyle', 'docked');
        subplot(6,1,1);
        PlotPsd(x, fs);
        subplot(6,1,2);
        PlotPsd(y, fs);
        subplot(6,1,3);
        PlotPsd(z, fs);
        subplot(6,1,4);
        hold on;
        stem(m, abs(zz), '.b');
        plot(m(i)*[1 1], get(gca, 'YLim'), '--r');
        subplot(6,1,5);
        hold on;
        stem(mm, abs(zzz), '.b');
        plot(mm(i)*[1 1], get(gca, 'YLim'), '--g');
        subplot(6,1,6);
        [zzzz, mmm] = xcorr(x, z, 100);
        stem(mmm, abs(zzzz), '.b');
    end
    
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
