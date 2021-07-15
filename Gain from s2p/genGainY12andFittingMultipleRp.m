clear variables
close all
clc

%% RESONATOR DATA
Qs = 8279;
C0 = 50e-15;
fs = 105.05e6;
oms = 2*pi*fs;
fp = 120.91e6;
kt2 = pi/2*fs/fp*1/tan(pi/2*fs/fp);

%% LOAD
%Load
Cl = 1000*1e-15;                     %Capacitive load [F]
RpV = [50 250 500 5000]*1e3;        %Resistance to ground at the output [Ohm]
Rs = 0;                             %Series resistance

%Source
Rin = 50;
Vs = 1;


%% Gain from Y12
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

for i = 1:length(RpV)

    Rp = RpV(i);
    %Calculates the gain provided by a series configuration
    Y = -y12;
    Yl = 1i*om*Cl+1/Rp;
    Zl = 1./Yl;

    Gs(i,:) = abs(Zl./(Rin+Rs+1./Y+Zl));

end

RpVs = string(RpV/1e6);
labelLoad = strcat("{\itR}_{p} = ",RpVs," M{\Omega}");

figure(1)
plot(freq/1e6, Gs,'LineWidth',3)
xlabel('Frequency, {\itf} [MHz]')
ylabel('Gain, {\itG} [V/V]')
legend(labelLoad)
set(gcf,'color','white')
set(gca,'FontSize',13)
grid on

axis([105.6 105.9 0 25])

%% Gain from fitting
% 
% %Frequency vector
% %freq = 0.8*fs:fs/10000:1.4*fs;
% omV = 2*pi*freq;
% 
% Rm = pi^2/8*1/(oms*Qs*kt2*C0); 
% Lm = pi^2/8*1/(oms^2*kt2*C0);
% Cm = 8/(pi^2)*kt2*C0;
% 
% Yr = 1i*omV*C0 + 1./(Rm + 1i*omV*Lm + 1./(1i*omV*Cm)); 
% 
% for i = 1:length(ClV)
% 
%     Cl = ClV(i);
%     %Gain
%     Zl = 1./(1i*omV*Cl);
% 
%     G(i,:) = abs(Zl./(Zl+1./Yr+Rin));
%     
% end
% 
% figure(2)
% plot(freq/1e6, G,'LineWidth',3)
% ylabel('Gain, {\itG} [V/V]')
% xlabel('Frequency, {\itf} [MHz]')
% legend(labelLoad)
% set(gcf,'color','w')
% set(gca,'FontSize',13)
% grid on
% 
% labelLoadM = strcat(labelLoad," (Meas.)");
% labelLoadF = strcat(labelLoad," (Fit.)");
% 
% %% Final plot
% figure(3)
% plot(freq/1e6, Gs,'LineWidth',3)
% hold on
% plot(freq/1e6, G,'-.','LineWidth',3)
% 
% xlabel('Frequency, {\itf} [MHz]')
% ylabel('Gain, {\itG} [V/V]')
% legend([labelLoadM labelLoadF])
% set(gcf,'color','white')
% set(gca,'FontSize',13)
% grid on