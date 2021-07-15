%clear all
%close all
%clc

%Load data
Vs = 1;
Rin = 50;
Cl = 1000e-15;

%Resonators data
fs = 50.84325e6;
oms = 2*pi*fs;
Q = 4064;
kt2 = 0.3094;
C0 = 113.47e-15;

%Rm = 26;
%Cm = 34.88e-15;
%Lm = 281.03e-6;

Rm = pi^2/8*1/(oms*Q*kt2*C0); 
Lm = pi^2/8*1/(oms^2*kt2*C0);
Cm = 8/(pi^2)*kt2*C0;


%Frequency vector
span_sx = 10e6;
span_dx = 20e6;
res = 1e3;

freqV = fs-span_sx:res:fs+span_dx;
omV = 2*pi*freqV;

Yr = 1i*omV*C0 + 1./(Rm + 1i*omV*Lm + 1./(1i*omV*Cm)); 

%Gain
Zl = 1./(1i*omV*Cl);

G = abs(Zl./(Zl+1./Yr+Rin));

figure(1)
plot(freqV/1e6, G,'LineWidth',3)
ylabel('Gain, G [V/V]')
xlabel('Frequency, f [MHz]')

set(gcf,'color','w')
set(gca,'FontSize',15)
grid on
