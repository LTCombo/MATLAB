function [dC0ei, dC0pp, C0ei, C0pp] = calc_C0_LFE(lambda, h, c, Le_lambda, Np, epsR)
    %Calculates the capacity density normalized over pair and electrode
    %length (dC0) with the elliptical integral (ei) and parallel plate (pp)
    %model 

    %Constants
    eps0 = 8.854187e-12;
    
    %% Elliptical Integral (EI)
    p = lambda./2;
    s = p*c;
    g = 0.5*p*(1-c);
    k0 = g/(g+s);
    k0p = sqrt(1-k0^2);
    k1 = tanh((pi*g)./(2*h))./tanh(pi*(s+g)./(2*h));
    k1p = sqrt(1-k1.^2);

    q1 = 0.5*(ellipke(k1p).*ellipke(k0))./(ellipke(k1).*ellipke(k0p));

    epseff = 1+(epsR-1).*q1;

    dC0ei = eps0.*epseff.*(ellipke(k0p)./ellipke(k0));
    
    %% Parallel plate model (PP)
    
    dC0pp = 2*eps0*epsR*h./(lambda.*(1-c));
    
    %% Capacitance calculation
    Le = Le_lambda*lambda;
    
    C0ei = dC0ei.*Le.*(Np+1);
    
    C0pp = dC0pp.*Le.*(Np+1);
    
end