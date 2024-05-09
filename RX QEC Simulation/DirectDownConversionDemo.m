 function [I_signal, Q_signal] = DirectDownConversionDemo()
    % Define parameters
    N = 1e5; % number of samples
    n = (0:N-1)'; % sample indices
    f_c = 0.1; % normalized carrier frequency
    phi_c = 0; % carrier phase
    phi_q = 0.0; % quadrature phase error
    g_q = 0.1; % quadrature gain imbalance

    % Generate a complex bandpass signal
    f_x = 0.05; % signal frequency
    x = exp(1i*2*pi*f_x*n); % complex exponential at f_x
    % % 
    % % Bandpass signal at RF with quadrature error
    % y = real(x .* exp(1i*(2*pi*f_c*n + phi_c)));
    % 
    % % Down-convert to baseband
    % y_i = y .* cos(2*pi*f_c*n + phi_c); % In-phase component
    % y_q = y .* sin(2*pi*f_c*n + phi_c); % Quadrature component
    % 
    % % Apply quadrature error
    g1 = (1 + g_q * exp(-1i*phi_q)) / 2;
    g2 = (1 - g_q * exp(1i*phi_q)) / 2;
    y_bb = x * exp(-1i*phi_c) .* g1 + conj(x) * exp(-1i*phi_c) .* g2; % Baseband signal with quadrature error
    
    % Apply low-pass filtering
    fbw = 0.2; % bandwidth of the LPF
    % y_bb_lpf = ApplyLpf(y_bb, fbw, true);

    % Sample the filtered signal (ADC)
    y_adc = y_bb; % 
    I_signal = real(y_adc);
    Q_signal = imag(y_adc);
    % Plot the spectrum of the baseband signal
    plot_ft(y_adc, 'Baseband Signal Spectrum with Quadrature Error');
end

function y = ApplyLpf(x, fbw, doPlot)
    
    % build filter
    M = 100;
    n = ((-M):(M))';
    w0 = 2*pi*fbw/2;
    h = sin(w0*n)./(pi*n);
    h(n == 0) = w0/pi;
    h = h .* kaiser(length(h), 10);
    
    % apply filter
    y = conv(x, h, 'valid');
    
    % We are simulating an analog filter, so there will still be
    % noise sources before we sample at the ADC
    % Let's add some noise back in to make the plots look nice.
    % (otherwise, the noise floor would get shaped by the filter too)
    y = y + 10^(-40/20)*randn(length(y), 1);
    
    % plot
    if ( (nargin >= 3) && doPlot )
        N = length(x);
        f = (-1/2):(1/N):(1/2-1/N);        
        wx = blackman(length(x));
        wy = blackman(length(y));
        X = fftshift(fft(wx.*x, N)) / sum(wx);
        Y = fftshift(fft(wy.*y, N)) / sum(wy);
        H = fftshift(fft(h, N));        
        figure;        
        subplot(2,1,1);
        hold on;
        plot(f, 20*log10(abs(X)), 'b');
        plot(f, 20*log10(abs(H)), 'r');
        legend('X(f)', 'H(f)');
        title('Input X(f)');
        xlabel('Frequency (normalized)');
        ylabel('Amplitude (dB)');
        xlim([-0.5 0.5]);
        ylim([-150 0]);
        subplot(2,1,2);
        hold on;
        plot(f, 20*log10(abs(Y)), 'b');
        plot(f, 20*log10(abs(H)), 'r');
        legend('Y(f)', 'H(f)');
        title('Output Y(f)');
        xlabel('Frequency (normalized)');
        ylabel('Amplitude (dB)');
        xlim([-0.5 0.5]);
        ylim([-150 0]);
    end
    
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
