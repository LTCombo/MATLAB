function Rs = calc_Rs_LFE(rho, tm, Lpad, Wpad, Lroute, La, Wa, W, bus, gap, c, Nf, Np, lambda, Le)

    p = lambda/2;

    Rpad = rho*Lpad/(2*Wpad*tm);
    Rroute = 2*rho*Lroute/((Wpad + Wa)*tm);
    Ra = rho*La/(Wa*tm);
    Rb = rho*W/(bus*tm);
    %Rf = 2*p*gap/(lambda*c*tm);
    Rf = rho*gap/(p*c*tm);
    Rg = 2*Rf/Nf;
    Re = 4/3*rho*Le/(lambda*c*tm)*1/Np;
    Rs = 2*(Rpad + Rroute + Ra + Rb + Rg) + Re;

end