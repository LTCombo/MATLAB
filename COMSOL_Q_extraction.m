% Preconditioning
clear
clc
close all

%% Load wide-span COMSOL simulation file

[y11, freq] = COMSOLreader();

y11_dB = 20*log10(abs(y11));
df = freq(2) - freq(1);

[y11Max, y11PosMax] = max(y11_dB);

y11Min3dB = y11Max-3;

BW3dB = (length( y11_dB(y11_dB > y11Min3dB)) + 1)*df;

Q3dB = freq(y11PosMax)/BW3dB;

figure(1)
plot(freq/1e6, 20*log10(abs(y11)))
xlabel('Frequency [MHz]')
ylabel('Y_{11} [dB]')

