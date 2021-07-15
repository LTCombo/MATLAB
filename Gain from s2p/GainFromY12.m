%% GAIN SIMULATOR FROM Y12
% By Luca Colombo - 02/13/2018
% lcolombo@andrew.cmu.edu
% Carnegie Mellon University

%Preconditioning
% clear variables
% close all
% clc
% 

%Source
Rin = 10;
Vs = 1;

%% Single run

%Load
Cl = 1000e-15;                   %Capacitive load [F]
Rp = 2e6;                     %Resistance to ground at the output [Ohm]
Rs = 0;                         %Series resistance

%Select the .s2p file for the Y12 extraction
[filename1, pathname] = uigetfile('Select file.s2p','Select the S file of series resonator');
filename = strcat(pathname, filename1);
data = read(rfdata.data, filename);                     %Note: requires MATLAB RF module

%Extracts freq and S parameters
freq = data.freq;
om=(2*pi).*freq;

s_params = extract(data, 'S_PARAMETERS',50);
s11 = squeeze(s_params(1,1,:));
s12 = squeeze(s_params(1,2,:));
s21 = squeeze(s_params(2,1,:));
s22 = squeeze(s_params(2,2,:));

%Converts to Y parameters
y_params = s2y(s_params, 50);
y11 = squeeze(y_params(1,1,:));
y12 = squeeze(y_params(1,2,:));
y21 = squeeze(y_params(2,1,:));
y22 = squeeze(y_params(2,2,:));

%Calculates the gain provided by a series configuration
Y = -y12;
Yl = 1i*om*Cl+1/Rp;
Zl = 1./Yl;

Gs = abs(Zl./(Rin+Rs+1./Y+Zl));

figure(1)
subplot(2,1,1)
plot(freq/1e6, Gs,'LineWidth',3)
xlabel('Frequency, {\itf} [MHz]')
ylabel('Gain, {\itG} [V/V]')

set(gcf,'color','white')
set(gca,'FontSize',15)
grid on

subplot(2,1,2)
plot(freq/1e6, 20*log10(Gs),'LineWidth',3)
xlabel('Frequency, {\itf} [MHz]')
ylabel('Gain, {\itG} [dB]')

set(gcf,'color','white')
set(gca,'FontSize',15)
grid on

%% Span

Cl_v = 10e-15:10e-15:1000e-15;                   %Capacitive load [F]
Rp = 2e6;                     %Resistance to ground at the output [Ohm]
Rs = 0;                         %Series resistance

for i = 1:length(Cl_v)
   
    Cl = Cl_v(i);
    Yl = 1i*om*Cl+1/Rp;
    Zl = 1./Yl;

    Gs = abs(Zl./(Rin+Rs+1./Y+Zl));
    [GsMax(i), posGsMax(i)] = max(Gs); 
    
end

figure(2)
subplot(2,1,1)
plot(Cl_v/1e-12, 20*log10(GsMax),'LineWidth',3)
xlabel('Capacitive load, {\itC}_{load} [pF]')
ylabel('Gain, {\itG} [dB]')

set(gcf,'color','white')
set(gca,'FontSize',15)
grid on

subplot(2,1,2)
plot(Cl_v/1e-12,freq(posGsMax)/1e6,'LineWidth',3)
xlabel('Capacitive load, {\itC}_{load} [pF]')
ylabel('Frequency, {\itf} [MHz]')

set(gcf,'color','white')
set(gca,'FontSize',15)
grid on
