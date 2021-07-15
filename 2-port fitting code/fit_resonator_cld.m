function [Q_3dB, Q_slope, kt2_b, fs, Rm, Lm, Cm, C0] = fit_resonator_cld(pathname, filename, name_res)

data = read(rfdata.data, char(filename));

%% read freq and s-parameters

    freq = data.freq;
    w=(2*pi).*freq;

    s_params = extract(data, 'S_PARAMETERS',50);
    s11 = squeeze(s_params(1,1,:));
    s12 = squeeze(s_params(1,2,:));
    s21 = squeeze(s_params(2,1,:));
    s22 = squeeze(s_params(2,2,:));

    % for reading y11, y12, y21, y22 from y-parameter
    y_params = s2y(s_params, 50);
    y11 = squeeze(y_params(1,1,:));
    y12 = squeeze(y_params(1,2,:));
    y21 = squeeze(y_params(2,1,:));
    y22 = squeeze(y_params(2,2,:));


%% converting Y to Z
z11 = 1./y11;

%impedance_plot=figure(1);
%semilogy(freq,(real(z11)), 'r'), 
%title('Real Impedance'), legend('Z11');

%% Plot of all frequency range 50MHz to 3GHz
%figure;
%plot_admittance_AlScN(freq, y11, data, 'r');

%% Fit Y11
FigName = strcat(strtok(name_res,'.'),' - Y11 Admittance Response without De-embedding');

%%
% figure;
% plot(f,imag(Y));
% figure;
% plot(f,real(Y));
samples = 100;
visible = 1 ;
have_frames = 0;

%% Some parameter that I may need for the extraction
f = freq;
delta_f = f(2) - f(1);
f_MHz = f*1e-6;
w = 2*pi*f;
Y = -y12;
Z = 1./Y;
YdB = 20*log10(abs(Y));
Yreal = real(Y);
Yimag = imag(Y);
len = length(Y); % this is the sample size
%samples = 2000; % this is the number of samples used for C0 and R0
                 % extraction; I can change this number accordingly to the 
                 % sample size

               
%% Rs
% This is the matching impedence of the resonator and is set to 50 Ohm.
Rs = 0;


%% C0
% This is the static capacitance of the resonator. I are going to
% calculate it from the out of band admitttance. I are going to take
% the avarage value obtained from the left and the right tail of the curve.
% I are going to use the immaginary part for the admittance.
% Mainly I use the definition X_C = 1/imag(Y) = 1/(j*w*C).
C01 = (mean(Yimag(1:samples)./w(1:samples))); % left tail
C02 = (mean(Yimag(len-samples:len)./w(len-samples:len))); % right tail
C0 = (C01+C02)/2; % medium value


%% R0
% This is the dielectric resistance of the resonator. I are going to
% calculate it from the out of band admittance. I are going to take
% the avarage value obtained from the left and the right tail of the curve.
% I are going to use the real part of the admittance.
% Mainly I use the definition X_R = 1/real(Y) = R.
Zreal = real(Z);
R1 = (mean(Zreal(1:samples))); % left tail
R2 = (mean(Zreal(end-samples:end))); % right tail
R = (R1+R2)/2; % medium value
R0 = R - Rs; % the two resistances are in serie


%% fs
% Here I are going to extract the resonance frequency from the serie
% impedence Zs. Basically I have the following equations:
% -> Z = Rs + Zs
% -> Z0 = R0 - 1/(j*w*C0)
% -> Zs = Zm // Z0 <==> 1/Zs = 1/Zm + 1/Z0 <==> Zm = (1/Zs - 1/Z0)^-1
% The resonance will be located at the minimum magnitude value of the 
% impedence Zm.
Zs = Z - Rs;
Z0 = R - 1./(1i*w*C0);
Zm = (1./Zs - 1./Z0).^-1;
[Zm_min_mag,fs_index] = min(abs(Zm));
fs = f(fs_index);
% This same result can be obtain directly from the resonator admittance Y
% and bypassing those calculations. The resonance will be located at the
% maximum magnitude value of the admittance Y. I do this just for check.
% The result will be slightly shifted so the first method is more accurate.
[Y_max,fs_index_check] = max(abs(Y));
fs_check = f(fs_index_check);

fs = fs_check;
fs_index = fs_index_check;


%% Rm
% This is the motional resistance. I extract it from the real part of
% the resonator admittance Y at his maximum amplitude or rather at
% resonance. In fact, at resonance, the capacitance in parallel with the
% inductance can be seen as a short circuit, so the only variable that
% matters is the serie resistance Rm.
Y_max_real = Yreal(fs_index);
R = 1 / Y_max_real;
Rm = R - Rs;
% Now we can do the same but use just Zm.
Zm_min_real = real(Zm(fs_index));
Rm_check = Zm_min_real;

Rm = max(1/abs(Y)) - Rs;

%% Q_3dB
% Here I'm going to calculate the Q based on the -3dB definition.
% The folmula would be: Q_3dB = fs / BWD_3dB.
Y_min_dB = min(YdB);
Y_max_dB = max(YdB);
Y_max_3dB = Y_max_dB - 3;
BWD = (length(YdB(YdB>Y_max_3dB))+1)*delta_f;
Q_3dB = fs/BWD;


%% Q_slope
% Here the Q is going to be computed as a frequency-to-phase slope.
n = 0;
h = f(fs_index+n) - f(fs_index+n-1);
ih = 1/(11*h);
df_dY = -180/pi*ih*(angle(Y(fs_index+n-2))-8*angle(Y(fs_index+n-1))+8*angle(Y(fs_index+n+1))-angle(Y(fs_index+n+2)));
y_intercept = 180/pi*angle(Y(fs_index+n)) + df_dY*f(fs_index+n);
% y = df_dY.*f + y_intercept;
Q_slope = pi*fs*df_dY/360;


%% Cm, Lm, Qm, fs
% This parameters can be computed from the resonator MBVD model circuit by 
% using the resonance's properties. We have the following equations:
[Y_min, fp_index] = min(abs(Y));
fp = f(fp_index);
% Backward calculations of some parameters.
Cm = C0*((fp/fs)^2-1);
Lm = 1/(Cm*(2*pi*fs)^2);
Q_m_recalc = 1/Rm*sqrt(Lm/Cm);


%% Some kt2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kt2_a = pi^2/8*Cm/C0;
kt2_b = pi^2/4*(1-fs/fp);
kt2_c = pi^2/8*(1-(fs/fp)^2);


%% Model fitting 1 @ no error correction
R0_est = 0;
Yfit_p = (1./(Rm+1i*w*Lm-1i./(w*Cm)))+(1./(1i*w*C0)+R0_est).^-1;
Zfit_p = 1./Yfit_p + Rs;
Yfit = 1./Zfit_p;
Yfit_dB = 20*log10(abs(Yfit));


%% Error detection
Err = Yfit - Y;
ErrImag = imag(Err);
ErrReal = real(Err);

%% Error on C0
Err_C01 = (mean(ErrImag(1:samples)./w(1:samples))); % left tail
Err_C02 = (mean(ErrImag(len-samples:len)./w(len-samples:len))); % right tail
Err_C0 = (Err_C01+Err_C02)/2; % medium value

%% Error on R0
Err_R1 = (mean(ErrReal(1:samples))); % left tail
Err_R2 = (mean(ErrReal(end-samples:end))); % right tail
Err_R = (Err_R1+Err_R2)/2; % medium value
Err_R0 = Err_R - Rs; % the two resistances are in serie


%% Model fitting 2 @ Err_C0
R0_est = 0;
C0 = C0 - Err_C0;

Cm = C0*((fp/fs)^2-1);
Lm = 1/(Cm*(2*pi*fs)^2);
Q_m_recalc_2 = 1/Rm*sqrt(Lm/Cm);
Q_m_recalc_2_Rs = 1/(Rm+Rs)*sqrt(Lm/Cm);

Yfit_p_2 = (1./(Rm+1i*w*Lm-1i./(w*Cm)))+(1./(1i*w*C0)+R0_est).^-1;
Zfit_p_2 = 1./Yfit_p_2 + Rs;
Yfit_2 = 1./Zfit_p_2;
Yfit_dB_2 = 20*log10(abs(Yfit_2));


%% Model fitting 3 @ Err_R0
R0 = R0 - Err_R;

Yfit_p_3 = (1./(Rm+1i*w*Lm-1i./(w*Cm)))+(1./(1i*w*C0)+R0).^-1;
Zfit_p_3 = 1./Yfit_p_3 + Rs;
Yfit_3 = 1./Zfit_p_3;
Yfit_dB_3 = 20*log10(abs(Yfit_3));

    
%% Plotting results
display_gap = (YdB(fs_index) - YdB(fp_index))*0.15;

%% New figure
if(visible)
    admittance_plot = figure;
else
    admittance_plot = figure('visible','off');
end

%% Admittance magnitude
if(have_frames)
    subplot(2,2,1);
else
    subplot(2,1,1);
end

plot(f_MHz, YdB, 'linewidth', 3); hold on;
%plot(f_MHz,Yfit_dB);
plot(f_MHz,Yfit_dB_2, 'linewidth', 3);
%plot(f_MHz,Yfit_dB_3);

Yfit_min_dB_2 = min(Yfit_dB_2);

xlim([f_MHz(1) f_MHz(end)]);
ylim([Yfit_min_dB_2-2 Y_max_dB+2]);

%title(FigName,'FontSize',16);

legend('Measured','Fitting','Location','SouthWest');

xlabel('Frequency, {\itf} [MHz]','FontSize',16);
ylabel('Amplitude, {\itY}_{12} [dB]','FontSize',16);

set(gca,'FontSize',16)
set(gcf,'color','white')
grid on



text(f_MHz(1)+(f_MHz(end)-f_MHz(1))*0.05, Y_max_dB-1*display_gap, strcat('f_s=', num2str(round(fs*1e-6,2)), ' MHz'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(1)+(f_MHz(end)-f_MHz(1))*0.05, Y_max_dB-2*display_gap, strcat('f_p=', num2str(round(fp*1e-6,2)), ' MHz'), 'FontSize',14,'FontWeight','Normal');

%text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-1*display_gap, strcat('R_s=', num2str(round(Rs,2)), ' \Omega'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-1*display_gap, strcat('C_0=', num2str(round(C0*1e15,2)), ' fF'), 'FontSize',14,'FontWeight','Normal');
%text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-3*display_gap, strcat('R_0=', num2str(round(R0_est,2)), ' \Omega'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-2*display_gap, strcat('R_m=', num2str(round(Rm,2)), ' \Omega'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-3*display_gap, strcat('L_m=', num2str(round(Lm*1e6,2)), ' uH'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-4*display_gap, strcat('C_m=', num2str(round(Cm*1e15,2)), ' fF'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-5*display_gap, strcat('k_t^2=', num2str(round(kt2_b*100,2)), ' %'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, Y_max_dB-6*display_gap, strcat('Q_{slope}=', num2str(round(Q_slope,2)), ''), 'FontSize',14,'FontWeight','Normal');

%% Admittance phase
if(have_frames)
    subplot(2,2,3);
else
    subplot(2,1,2);
end

plot(f_MHz,angle(Y)/pi*180, 'linewidth', 3); hold on;
plot(f_MHz,angle(Yfit_2)/pi*180, 'linewidth', 3);

legend('Measured','Fitting','Location','SouthWest');

xlim([f_MHz(1) f_MHz(end)]);
ylim([-100 100]);

xlabel('Frequency, {\itf} [MHz]','FontSize',16);
ylabel('Phase, \phi [\circ]','FontSize',16);

set(gca,'FontSize',16)
set(gcf,'color','white')
grid on

x0=100;
y0=100;

width=1000;
height=800;
set(gcf,'position',[x0,y0,width,height])

text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, 80-2*display_gap, strcat('k_t^2=', num2str(round(kt2_b*100,2)), ' %'), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, 80-5*display_gap, strcat('Q_{slope}=', num2str(round(Q_slope,2)), ''), 'FontSize',14,'FontWeight','Normal');
%text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, 80-8*display_gap, strcat('Q_{m_{Rs}}=', num2str(round(Q_m_recalc_2_Rs,2)), ''), 'FontSize',14,'FontWeight','Normal');
text(f_MHz(end)-(f_MHz(end)-f_MHz(1))*0.25, 80-8*display_gap, strcat('Qm_{Err_{C0}}=', num2str(round(Q_m_recalc_2,2)), ''), 'FontSize',14,'FontWeight','Normal');

%% Plot legend
%hold off;
%legend('Measured','Fitting','Location','SouthWest');
%% Grid on
grid on;

%%Save figures

label_1='_impedance';
label_2='_admittance';
ext_1='.jpg';
ext_2='.fig';

name_res=erase(name_res,'.s2p');
admittance_name_jpg=strcat(strcat(pathname,'Figures\'),strcat(char(name_res),ext_1));
admittance_name_fig=strcat(strcat(pathname,'Figures\'),strcat(char(name_res),ext_2));

saveas(admittance_plot,admittance_name_jpg);
saveas(admittance_plot,admittance_name_fig);

close all

end