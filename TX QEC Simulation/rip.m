function P = rip(gq1, tq1, gq2, tq2)

    % assume gq1,tq1 represent the actual analog impairments
    g1 = 0.5*(1 + gq1.*exp(1i*tq1));
    g2 = 0.5*(1 - gq1.*exp(1i*tq1));

    % assume gq2,tq2 are used for correction
    b = 1/(gq2*cos(tq2)) - 1i*tan(tq2);
    
    % compute residual image power
    P = abs( (g1.*(1-b) + g2.*(1+conj(b))) ./ (g1.*(1+b) + g2.*(1-conj(b))) ).^2;
    
end
