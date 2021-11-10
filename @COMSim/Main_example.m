
%%Preconditioning
clear batch
close all
clc

%Help
%mphdoc(batch.model.study)

%batch=preamble();

%% Data
freq = 10:10:300;
study = 'std5';

%% Model Inizialiation
batch = COMSim();
batch.study = batch.model.study(study);
batch.set_freqs_with_unit(study, freq, 'MHz');

% hol Sweep
batch.sweep_param='hol';
%batch.sweep_unit='[]';
batch.sweep_value=0.01:0.01:0.5;

batch.setup_folder();
batch.run_loop();
%batch.analyze_data();
%batch.delete;

%% Plot
figure(1)
plot(freq, 20*log10(abs(batch.data')), 'LineWidth', 3)
xlabel('Frequency, {\itf} [MHz]')
ylabel('Admittance, {\itY} [dB]')

plotPreferences

%% Save
data.freq = freq;
data.param = batch.sweep_param;
data.param_sweep = batch.sweep_value;
data.Y = batch.data;

save("hol_C0_sweep.mat", 'data')
