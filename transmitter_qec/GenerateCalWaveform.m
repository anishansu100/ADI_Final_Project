function x = GenerateCalWaveform(filename, gq, tq)

    if (nargin < 2)
        gq = 1;
        tq = 0;
    end

    % number of samples to generate
    % PLUTO Tx/Rx buffers are 64k samples long
    N = 2^16;
    
    % sampling rate in Hz
    FS = 30.72e6;
    
    % sample indices
    n = (0:(N-1))';
    
    % bandwidth
    B = 500e3;
    
    % frequency offset
    C = 2.5e6;
    
    % PAPR target
    PAPR_TARGET = 10^(4/10);
    
    % generate a narrowband modulated waveform    
    U = ceil(FS/B);
    M = ceil(N/U);
    x = randn(M,1) + 1i*randn(M,1);
    [x, b] = resample(x, U, 1);
    x = x(1:N) .* exp(1i*2*pi*C/FS*n);    
    b = b(:) .* exp(1i*2*pi*C/FS*(1:length(b))');
    
    % apply CFR
    for i = 1:4
        x = RunCfr(x, PAPR_TARGET, b);
    end
    
    % scale
    x = 10^(-1/20) * x / max(max(abs(real(x))), max(abs(imag(x))));
    
    % apply QEC correction
    x = (real(x) + tan(tq)*imag(x)) + 1i*(imag(x)/(gq*cos(tq)));
    
    % quantize to 12 bits
    x = round(x * 2^11) / 2^11;        
    
    % write waveform to file that can be read by Python and loaded into Tx
    if (nargin > 0 && ~isempty(filename))
        csvwrite(filename, [real(x(:)) imag(x(:))]);
    end
    
    if (0) % debug
        figure;
        set(gcf, 'WindowStyle', 'docked');    
        subplot(4,1,1);
        hold on;
        plot(real(x), 'b');
        plot(imag(x), 'r');
        title(sprintf('RMS = %f dBFS, PAPR = %f dB', 10*log10(mean(x.*conj(x))), 10*log10(GetPAPR(x))));
        subplot(4,1,2);
        plot(abs(x), 'b');
        subplot(4,1,3);
        PlotPsd(x, FS);
        subplot(4,1,4);    
        [xx, m] = xcorr(x, 100);
        stem(m, abs(xx), '.b');
    end

end

function papr = GetPAPR(x)
    p = max(x.*conj(x));
    a = mean(x.*conj(x));
    papr = p/a;
end

function y = RunCfr(x, PAPR_TARGET, b)
    A = sqrt(mean(x.*conj(x)) * PAPR_TARGET);
    if (0)
        % clipping
        y = x;
        m = abs(y) > A;
        y(m) = (A./abs(y(m))) .* y(m);
    else
        % peak cancellation
        C = 0.9*A;
        y = x;
        if (mod(length(b),2) == 1)
            M = (length(b)-1)/2;
            n = (-M):M;
        else
            M = length(b)/2;
            n = (-M):(M-1);
        end
        while (1)
            [a,i] = max(abs(y));
            if (a < A)
                break;
            end            
            j = i + n;
            m = and(j >= 1, j <= length(y));
            
            % Amplitude correction
            % |y - c.b|^2 = C^2
            % (y - c.b)(y* - c*.b*) = C^2
            % |y|^2 - y.c*.b* - y*.c.b + |c|^2.|b|^2 = C^2            
            % |y|^2 - y.(ci-j.cq).b* - y*.(ci+j.cq).b + (ci^2+cq^2).|b|^2 = C^2
            % (-y.b*-y*.b).ci + j.(y.b*-y*.b).cq + (|b|^2).ci^2 + (|b|^2).cq^2 = C^2 - |y|^2
            % (-2Re{yb*}).ci + (-2Im{yb*}).cq + (|b|^2).ci^2 + (|b|^2).cq^2 = C^2 - |y|^2
            % |b|^2.ci^2 + |b|^2.cq^2 - 2Re{yb*}.ci - 2Im{yb*}.cq = C^2 - |y|^2
            %
            % Phase correction (keep same phase)
            % ci/cq = yi/yq ==> cq = ci*yq/yi            
            %
            % Put it together
            % [|b|^2(1+yq^2/yi^2)].ci^2 + [-(2Re{yb*}+2Im{yb*}.yq/yi)].ci + [|y|^2 - C^2] = 0
            % D2.ci^2 + D1.ci + D0 = 0
            D2 = abs(b(n==0))^2*(1 + imag(y(i))^2/real(y(i))^2);
            D1 = -1*(2*real(y(i)*conj(b(n==0))) + 2*imag(y(i)*conj(b(n==0)))*imag(y(i))/real(y(i)));
            D0 = abs(y(i))^2 - C^2;
            ci = (-D1 + [-1 1]*sqrt(D1^2 - 4*D2*D0)) / (2*D2);
            cq = ci*imag(y(i))/real(y(i));
            [~,k] = min(abs(ci+1i*cq));
            c = ci(k) + 1i*cq(k);            
            
            y(j(m)) = y(j(m)) - c*b(m);            
        end
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
