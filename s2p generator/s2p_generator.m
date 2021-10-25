%MEMS Resonator s2p generator
%Luca Colombo (l.colombo@northeastern.edu) - 08/02/2021

%Preconditioning
clear all
close all
clc

%Design parameters
fs = 10e6;                  %Center frequency [Hz]
Qs = 10000;                 %Quality factor
kt2 = 0.3;                  %EM coupling
Cload = 200e-15;            %ED input capacitance [F]
Rload = 1e6;                %ED input resistance [Ohm]
C0 = Cload;                 %Resonator static capacitance [F]
Cp = 15e-15;                %Capacitive parasitic to ground [F] 
Rp = 1e6;                   %Resistive parasitic to ground [Ohm]
Cf = 5e-15;                 %Feedthrough capacitance [F]
Rf = 10e6;                  %Feedthrough parallel resistance [Ohm]

%Port and simulations parameters
Zin = 50;                                   %Source impedance [Ohm]
fStart = 0.7*fs;                            %Starting frequency [Hz]
fEnd = 1.45*fs;                              %Ending frequency [Hz]
nPoints = 1000;                             %Number of simulated points
freq = linspace(fStart, fEnd, nPoints);     %Frequency vector [Hz]
om = 2*pi*freq;                             %Angular frequency [rad/s]

%Resonator model
Rm = pi^2/8*1/(2*pi*fs*Qs*kt2*C0);          %Motional resistance [Ohm]
Lm = pi^2/8*1/((2*pi*fs)^2*kt2*C0);         %Motional inductance [H]
Cm = 8/(pi^2)*kt2*C0;                       %Motional capacitance [F]

%Building BVD model
%cfr. Piazza Piezoelectric Resonator book for more explanations

Y = 1i*om*C0 + 1./(Rm + 1i*om*Lm + 1./(1i*om*Cm));

figure(1)
plot(freq/1e6, 20*log10(abs(Y)), 'LineWidth', 3)
grid on
xlabel('Frequency, {\itf} [MHz]')
ylabel('Admittance, {\itY} [dB]')
set(gcf,'color','white')
set(gca,'fontSize',13)

%Adding parasitics
Yf = 1i*om*Cf + 1./(Rf);                    %Feedthrough admittance [S]
Yp = 1i*om*Cp + 1./(Rp);                    %Parasitic-to-ground admittance [S]

Y11 = Y + Yf + Yp; 
Y22 = Y11;
Y12 = -(Y + Yf);
Y21 = Y12;

%Tabulate Y-parameters
Y_param(1,1,:) = Y11;
Y_param(1,2,:) = Y12;
Y_param(2,1,:) = Y21;
Y_param(2,2,:) = Y22;

%Generate S-parameters
S_param = y2s(Y_param, Zin);

%Write Touchstone file
rfwrite(S_param, freq, 'data.s2p');