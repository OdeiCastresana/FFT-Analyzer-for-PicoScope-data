% Read the data from both CSV files, skipping the first 3 rows (header)
data_pm1 = dlmread('130520250_080PF_100pct_PM1_no_filt.csv', ',', 3, 0);
data_pm2 = dlmread('140520250_080PF_100pct_PM2_cap.csv', ',', 3, 0);

% Confirm data loaded
disp('Data for PM1 and PM2 loaded successfully');

% Separate columns for PM1
time_pm1 = data_pm1(:, 1);
voltage_pm1 = data_pm1(:, 3);
voltage_cap_pm1 = data_pm1(:, 2);

% Separate columns for PM2
time_pm2 = data_pm2(:, 1);
voltage_pm2 = data_pm2(:, 3);
voltage_cap_pm2 = data_pm2(:, 2);

% Resistances for PM1 and PM2 (ohms)
R_pm1 = 0.4181;
R_pm2 = 0.3573;

% Current for PM1 from Ohm's law
current_pm1 = voltage_pm1 / R_pm1;

% Current for PM2 from Ohm's law
current_pm2 = voltage_pm2 / R_pm2;

% Sampling frequency (given): 80 MHz
Fs = 80e6;
disp(['Sampling Frequency Fs = ', num2str(Fs), ' Hz']);

% Number of points for PM1
n_pm1 = length(voltage_pm1);

% Number of points for PM2
n_pm2 = length(voltage_pm2);

% FFT computation for PM1
fft_voltage_pm1 = fft(voltage_pm1);
fft_current_pm1 = fft(current_pm1);
fft_voltage_cap_pm1 = fft(voltage_cap_pm1);

% FFT computation for PM2
fft_voltage_pm2 = fft(voltage_pm2);
fft_current_pm2 = fft(current_pm2);
fft_voltage_cap_pm2 = fft(voltage_cap_pm2);

% Frequency vector for PM1
frequencies_pm1 = (0:n_pm1-1)*(Fs/n_pm1);

% Frequency vector for PM2
frequencies_pm2 = (0:n_pm2-1)*(Fs/n_pm2);

% Magnitudes normalized for PM1
mag_voltage_pm1 = abs(fft_voltage_pm1) * 2 / n_pm1;
mag_current_pm1 = abs(fft_current_pm1) * 2 / n_pm1;
mag_voltage_cap_pm1 = abs(fft_voltage_cap_pm1) * 2 / n_pm1;

% Magnitudes normalized for PM2
mag_voltage_pm2 = abs(fft_voltage_pm2) * 2 / n_pm2;
mag_current_pm2 = abs(fft_current_pm2) * 2 / n_pm2;
mag_voltage_cap_pm2 = abs(fft_voltage_cap_pm2) * 2 / n_pm2;

% Filter to 10 Hz – 10 MHz for both PM1 and PM2
low_lim = 10;
high_lim = 5e7;
idx_pm1 = (frequencies_pm1 >= low_lim) & (frequencies_pm1 <= high_lim);
freq_filt_pm1 = frequencies_pm1(idx_pm1);
mag_voltage_filt_pm1 = mag_voltage_pm1(idx_pm1);
mag_current_filt_pm1 = mag_current_pm1(idx_pm1);
mag_voltage_cap_filt_pm1 = mag_voltage_cap_pm1(idx_pm1);

idx_pm2 = (frequencies_pm2 >= low_lim) & (frequencies_pm2 <= high_lim);
freq_filt_pm2 = frequencies_pm2(idx_pm2);
mag_voltage_filt_pm2 = mag_voltage_pm2(idx_pm2);
mag_current_filt_pm2 = mag_current_pm2(idx_pm2);
mag_voltage_cap_filt_pm2 = mag_voltage_cap_pm2(idx_pm2);

% Plot PM1 (Voltage and Current FFT)
% Plot PM1 (Voltage, Current, Capacitor Voltage FFT)
figure;
subplot(3,1,1);
semilogx(freq_filt_pm1, mag_voltage_filt_pm1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Voltage Amplitude (V)');
title('PM1 - Snubber Resistor Voltage FFT');
xlim([low_lim high_lim]);
ylim([0 max(mag_voltage_filt_pm1) + 0.1]);

subplot(3,1,2);
semilogx(freq_filt_pm1, mag_current_filt_pm1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Current Amplitude (A)');
title('PM1 - Snubber Resistor Current FFT');
xlim([low_lim high_lim]);
ylim([0 max(mag_current_filt_pm1) + 0.1]);

subplot(3,1,3);
semilogx(freq_filt_pm1, mag_voltage_cap_filt_pm1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Voltage Amplitude (V)');
title('PM1 - Snubber Capacitor Voltage FFT');
xlim([low_lim high_lim]);
ylim([0 max(mag_voltage_cap_filt_pm1) + 0.1]);

print('PM1_FFT_Plots.png', '-dpng');  % Save PM1 figure as PNG

% Plot PM2 (Voltage, Current, Capacitor Voltage FFT)
figure;
subplot(3,1,1);
semilogx(freq_filt_pm2, mag_voltage_filt_pm2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Voltage Amplitude (V)');
title('PM2 - Snubber Resistor Voltage FFT');
xlim([low_lim high_lim]);
ylim([0 max(mag_voltage_filt_pm2) + 0.1]);

subplot(3,1,2);
semilogx(freq_filt_pm2, mag_current_filt_pm2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Current Amplitude (A)');
title('PM2 - Snubber Resistor Current FFT');
xlim([low_lim high_lim]);
ylim([0 max(mag_current_filt_pm2) + 0.1]);

subplot(3,1,3);
semilogx(freq_filt_pm2, mag_voltage_cap_filt_pm2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Voltage Amplitude (V)');
title('PM2 - Snubber Capacitor Voltage FFT');
xlim([low_lim high_lim]);
ylim([0 max(mag_voltage_cap_filt_pm2) + 0.1]);

print('PM2_FFT_Plots.png', '-dpng');  % Save PM2 figure as PNG


% RMS, RMS AC, DC, Max, Min for PM1
dc_v_pm1 = mean(voltage_pm1);
dc_i_pm1 = mean(current_pm1);
dc_v_cap_pm1 = mean(voltage_cap_pm1);
rms_v_pm1 = sqrt(mean(voltage_pm1.^2));
rms_i_pm1 = sqrt(mean(current_pm1.^2));
rms_v_cap_pm1 = sqrt(mean(voltage_cap_pm1.^2));
rms_ac_v_pm1 = sqrt(mean((voltage_pm1 - dc_v_pm1).^2));
rms_ac_i_pm1 = sqrt(mean((current_pm1 - dc_i_pm1).^2));
rms_ac_v_cap_pm1 = sqrt(mean((voltage_cap_pm1 - dc_v_cap_pm1).^2));
max_v_pm1 = max(voltage_pm1); min_v_pm1 = min(voltage_pm1);
max_i_pm1 = max(current_pm1); min_i_pm1 = min(current_pm1);
max_v_cap_pm1 = max(voltage_cap_pm1); min_v_cap_pm1 = min(voltage_cap_pm1);

% RMS, RMS AC, DC, Max, Min for PM2
dc_v_pm2 = mean(voltage_pm2);
dc_i_pm2 = mean(current_pm2);
dc_v_cap_pm2 = mean(voltage_cap_pm2);
rms_v_pm2 = sqrt(mean(voltage_pm2.^2));
rms_i_pm2 = sqrt(mean(current_pm2.^2));
rms_v_cap_pm2 = sqrt(mean(voltage_cap_pm2.^2));
rms_ac_v_pm2 = sqrt(mean((voltage_pm2 - dc_v_pm2).^2));
rms_ac_i_pm2 = sqrt(mean((current_pm2 - dc_i_pm2).^2));
rms_ac_v_cap_pm2 = sqrt(mean((voltage_cap_pm2 - dc_v_cap_pm2).^2));
max_v_pm2 = max(voltage_pm2); min_v_pm2 = min(voltage_pm2);
max_i_pm2 = max(current_pm2); min_i_pm2 = min(current_pm2);
max_v_cap_pm2 = max(voltage_cap_pm2); min_v_cap_pm2 = min(voltage_cap_pm2);

% Fundamental (peak in spectrum) for PM1
[~, idx_v_pm1] = max(mag_voltage_filt_pm1);
fund_v_pm1 = freq_filt_pm1(idx_v_pm1);
amp_v_pm1 = mag_voltage_filt_pm1(idx_v_pm1);
phase_v_pm1 = angle(fft_voltage_pm1(idx_v_pm1));

[~, idx_i_pm1] = max(mag_current_filt_pm1);
fund_i_pm1 = freq_filt_pm1(idx_i_pm1);
amp_i_pm1 = mag_current_filt_pm1(idx_i_pm1);
phase_i_pm1 = angle(fft_current_pm1(idx_i_pm1));

[~, idx_v_cap_pm1] = max(mag_voltage_cap_filt_pm1);
fund_v_cap_pm1 = freq_filt_pm1(idx_v_cap_pm1);
amp_v_cap_pm1 = mag_voltage_cap_filt_pm1(idx_v_cap_pm1);
phase_v_cap_pm1 = angle(fft_voltage_cap_pm1(idx_v_cap_pm1));

% Fundamental (peak in spectrum) for PM2
[~, idx_v_pm2] = max(mag_voltage_filt_pm2);
fund_v_pm2 = freq_filt_pm2(idx_v_pm2);
amp_v_pm2 = mag_voltage_filt_pm2(idx_v_pm2);
phase_v_pm2 = angle(fft_voltage_pm2(idx_v_pm2));

[~, idx_i_pm2] = max(mag_current_filt_pm2);
fund_i_pm2 = freq_filt_pm2(idx_i_pm2);
amp_i_pm2 = mag_current_filt_pm2(idx_i_pm2);
phase_i_pm2 = angle(fft_current_pm2(idx_i_pm2));

[~, idx_v_cap_pm2] = max(mag_voltage_cap_filt_pm2);
fund_v_cap_pm2 = freq_filt_pm2(idx_v_cap_pm2);
amp_v_cap_pm2 = mag_voltage_cap_filt_pm2(idx_v_cap_pm2);
phase_v_cap_pm2 = angle(fft_voltage_cap_pm2(idx_v_cap_pm2));

% Save results for PM1
output_file_pm1 = 'PM1_results.txt';
fid_pm1 = fopen(output_file_pm1, 'w');
fprintf(fid_pm1, 'Data for Capacitor Life Estimation Analysis for PM1 (1A33)\n');
fprintf(fid_pm1, '---------------------------------\n');

% Voltage (resistance) data for PM1
fprintf(fid_pm1, '\nSignal Voltage PM1: 1A33_V_R_snubber\n');
fprintf(fid_pm1, '---------------------------------\n');
fprintf(fid_pm1, 'DC: %.4f V\n', dc_v_pm1);
fprintf(fid_pm1, 'Max: %.4f V\n', max_v_pm1);
fprintf(fid_pm1, 'Min: %.4f V\n', min_v_pm1);
fprintf(fid_pm1, 'RMS: %.4f V\n', rms_v_pm1);
fprintf(fid_pm1, 'RMS AC: %.4f V\n', rms_ac_v_pm1);
fprintf(fid_pm1, 'Fundamental frequency: %.2f Hz\n', fund_v_pm1);
fprintf(fid_pm1, 'Amplitude is: %.4f V\n', amp_v_pm1);
fprintf(fid_pm1, 'Phase is: %.4f rad\n', phase_v_pm1);

% Separator for current
fprintf(fid_pm1, '\n--------------------------------------------------\n');
fprintf(fid_pm1, 'Signal Current PM1: 1A33_I_R_snubber\n');
fprintf(fid_pm1, '--------------------------------------\n');
fprintf(fid_pm1, 'DC: %.4f A\n', dc_i_pm1);
fprintf(fid_pm1, 'Max: %.4f A\n', max_i_pm1);
fprintf(fid_pm1, 'Min: %.4f A\n', min_i_pm1);
fprintf(fid_pm1, 'RMS: %.4f A\n', rms_i_pm1);
fprintf(fid_pm1, 'RMS AC: %.4f A\n', rms_ac_i_pm1);
fprintf(fid_pm1, 'Fundamental frequency: %.2f Hz\n', fund_i_pm1);
fprintf(fid_pm1, 'Amplitude is: %.4f A\n', amp_i_pm1);
fprintf(fid_pm1, 'Phase is: %.4f rad\n', phase_i_pm1);

% Separator for capacitor voltage
fprintf(fid_pm1, '\n--------------------------------------------------\n');
fprintf(fid_pm1, 'Signal Voltage PM1: 1A33_V_C_snubber\n');
fprintf(fid_pm1, '---------------------------------------\n');
fprintf(fid_pm1, 'DC: %.4f V\n', dc_v_cap_pm1);
fprintf(fid_pm1, 'Max: %.4f V\n', max_v_cap_pm1);
fprintf(fid_pm1, 'Min: %.4f V\n', min_v_cap_pm1);
fprintf(fid_pm1, 'RMS: %.4f V\n', rms_v_cap_pm1);
fprintf(fid_pm1, 'RMS AC: %.4f V\n', rms_ac_v_cap_pm1);
fprintf(fid_pm1, 'Fundamental frequency: %.2f Hz\n', fund_v_cap_pm1);
fprintf(fid_pm1, 'Amplitude is: %.4f V\n', amp_v_cap_pm1);
fprintf(fid_pm1, 'Phase is: %.4f rad\n\n', phase_v_cap_pm1);

% Save results for PM2
output_file_pm2 = 'PM2_results.txt';
fid_pm2 = fopen(output_file_pm2, 'w');
fprintf(fid_pm2, 'Data for Capacitor Life Estimation Analysis for PM2 (2A33)\n');
fprintf(fid_pm2, '---------------------------------\n');

% Voltage (resistance) data for PM2
fprintf(fid_pm2, '\nSignal Voltage PM2: 2A33_V_R_snubber\n');
fprintf(fid_pm2, '---------------------------------\n');
fprintf(fid_pm2, 'DC: %.4f V\n', dc_v_pm2);
fprintf(fid_pm2, 'Max: %.4f V\n', max_v_pm2);
fprintf(fid_pm2, 'Min: %.4f V\n', min_v_pm2);
fprintf(fid_pm2, 'RMS: %.4f V\n', rms_v_pm2);
fprintf(fid_pm2, 'RMS AC: %.4f V\n', rms_ac_v_pm2);
fprintf(fid_pm2, 'Fundamental frequency: %.2f Hz\n', fund_v_pm2);
fprintf(fid_pm2, 'Amplitude is: %.4f V\n', amp_v_pm2);
fprintf(fid_pm2, 'Phase is: %.4f rad\n', phase_v_pm2);

% Separator for current
fprintf(fid_pm2, '\n--------------------------------------------------\n');
fprintf(fid_pm2, 'Signal Current PM2: 2A33_I_R_snubber\n');
fprintf(fid_pm2, '--------------------------------------\n');
fprintf(fid_pm2, 'DC: %.4f A\n', dc_i_pm2);
fprintf(fid_pm2, 'Max: %.4f A\n', max_i_pm2);
fprintf(fid_pm2, 'Min: %.4f A\n', min_i_pm2);
fprintf(fid_pm2, 'RMS: %.4f A\n', rms_i_pm2);
fprintf(fid_pm2, 'RMS AC: %.4f A\n', rms_ac_i_pm2);
fprintf(fid_pm2, 'Fundamental frequency: %.2f Hz\n', fund_i_pm2);
fprintf(fid_pm2, 'Amplitude is: %.4f A\n', amp_i_pm2);
fprintf(fid_pm2, 'Phase is: %.4f rad\n', phase_i_pm2);

% Separator for capacitor voltage
fprintf(fid_pm2, '\n--------------------------------------------------\n');
fprintf(fid_pm2, 'Signal Voltage PM2: 2A33_V_C_snubber\n');
fprintf(fid_pm2, '---------------------------------------\n');
fprintf(fid_pm2, 'DC: %.4f V\n', dc_v_cap_pm2);
fprintf(fid_pm2, 'Max: %.4f V\n', max_v_cap_pm2);
fprintf(fid_pm2, 'Min: %.4f V\n', min_v_cap_pm2);
fprintf(fid_pm2, 'RMS: %.4f V\n', rms_v_cap_pm2);
fprintf(fid_pm2, 'RMS AC: %.4f V\n', rms_ac_v_cap_pm2);
fprintf(fid_pm2, 'Fundamental frequency: %.2f Hz\n', fund_v_cap_pm2);
fprintf(fid_pm2, 'Amplitude is: %.4f V\n', amp_v_cap_pm2);
fprintf(fid_pm2, 'Phase is: %.4f rad\n\n', phase_v_cap_pm2);

% Close the files
fclose(fid_pm1);
fclose(fid_pm2);

disp('Results saved to PM1_results.txt and PM2_results.txt');


