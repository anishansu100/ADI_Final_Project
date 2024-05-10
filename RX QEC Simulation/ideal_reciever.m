function ideal_reciever()
    % % number of samples to generate
    N = 1e5;
    % 
    % % sample indices
    n = (0:(N-1))';
    % 
    % % signal x(t)
    f_x = 0.05;
    % % Define parameters for the bandpass signal
    f_c = 0.1; % Normalized carrier frequency (as a fraction of the sampling rate)
    phi_c = 0; % Carrier phase
    % 
    % % generate a complex exponential waveform
    % % x[n] = e^(j*w0*n) = e^(j*2*pi*f0*n)    
    x = exp(1i*2*pi*f_x*n);
    
    % multiply the signal we want to decode with e^(j*2*pi*f_c*n)  
    y = x .* exp(1i*(2*pi*f_c*n + phi_c)); 
    y = real(y) + 10^(-50/20)*(randn(N,1));

    % % add some noise to make the plots look more realistic   
    % y = y + 10^(-50/20)*(randn(N,1)); 

    % Rx mixer down-converts to baseband
    y_i = y .* cos(2*pi*f_c*n + phi_c); % In-phase component
    y_q = -y .* sin(2*pi*f_c*n + phi_c); % Quadrature component
    
    % Combine in-phase and quadrature components to form complex baseband signal
    y_bb = y_i + 1i*y_q;
    
    % Apply low-pass filtering (simple moving average filter)
    % lpf_window_size = 100; % Window size of the moving average filter
    % y_bb_lpf = filter(ones(lpf_window_size,1)/lpf_window_size, 1, y_bb);
    y_bb_lpf = ApplyLpf(y_bb, 0.2, false);
    
    % FFT of the filtered signal
    Y_bb = y_bb_lpf;
    
   

    plot_ft(Y_bb, "Ideal Reciever with No Quadrature Error");
    % Plot the original and down-converted signal spectra
    % figure;
    % subplot(3,1,1);
    % plot(f, 20*log10(abs(fftshift(fft(x)))/N));
    % ylabel('Amplitude (dB)');
    % xlabel('Frequency (normalized)');
    % title('|X(f)|^2');
    % 
    % subplot(3,1,2);
    % plot(f, 20*log10(abs(fftshift(fft(y_bb)))/N));
    % ylabel('Amplitude (dB)');
    % xlabel('Frequency (normalized)');
    % title('|Y_{bb}(f)|^2 Before LPF');
    % 
    % subplot(3,1,3);
    % plot(f, 20*log10(abs(Y_bb)/N));
    % ylabel('Amplitude (dB)');
    % xlabel('Frequency (normalized)');
    % title('|Y_{bb}(f)|^2 After LPF');

end


function plot_ft(x, t)
    N = length(x);
    % Compute frequency vector for plotting
    f = (-N/2):(N-1)/2;
    f = f / N; % Normalize frequencies
    
    figure;
    plot(f, 20*log10(abs(fftshift(fft(x)))/N));
    ylabel('Amplitude (dB)');
    xlabel('Frequency (normalized)');
    title(t);

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
