% Read the data from the CSV file, skipping the first 3 rows (header)
data_pm1 = dlmread('130520250_08PF_80pct_PM2_cap.csv', ',', 3, 0);

% Confirm data loaded
disp('PM1 data loaded successfully');

% Separate columns for PM1
time_pm1 = data_pm1(:, 1);
voltage_pm1 = data_pm1(:, 3);
voltage_cap_pm1 = data_pm1(:, 2);

% Resistance for PM1 (ohms)
R_pm1 = 0.4181;

% Current for PM1 using Ohm's law
current_pm1 = voltage_pm1 / R_pm1;

% Sampling frequency (given): 80 MHz
Fs = 80e6;
disp(['Sampling frequency Fs = ', num2str(Fs), ' Hz']);

% Number of data points for PM1
n_pm1 = length(voltage_pm1);

% Compute FFT for PM1 voltage
fft_voltage_pm1 = fft(voltage_pm1);

% Frequency vector for PM1
frequencies_pm1 = (0:n_pm1-1) * (Fs / n_pm1);

% Normalize FFT magnitudes
mag_voltage_pm1 = abs(fft_voltage_pm1) * 2 / n_pm1;

% Filter frequencies between 10 Hz and 50 MHz
low_lim = 10;
high_lim = 5e7;
idx_pm1 = (frequencies_pm1 >= low_lim) & (frequencies_pm1 <= high_lim);

freq_filt_pm1 = frequencies_pm1(idx_pm1);
mag_filt_pm1 = mag_voltage_pm1(idx_pm1);

% Apply threshold: discard magnitudes lower than 0.1% of the maximum
max_mag_pm1 = max(mag_filt_pm1);
threshold = 0.001 * max_mag_pm1;  % 0.1% of max magnitude

idx_threshold = mag_filt_pm1 >= threshold;

freq_final = freq_filt_pm1(idx_threshold);
mag_final = mag_filt_pm1(idx_threshold);

% Export filtered FFT data to CSV
csvwrite('PM1_FFT_filtered.csv', [freq_final(:), mag_final(:)]);
disp('File PM1_FFT_filtered.csv created with filtered FFT above 0.1% threshold');

% === Plot the filtered FFT data ===
figure;
semilogx(freq_final, mag_final, 'b.-');
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (V)');
title('PM1 FFT - Components ≥ 0.1% of maximum amplitude');
xlim([low_lim high_lim]);
ylim([0 max(mag_final) * 1.1]);


