% E101: HW 9a
% Sebastian Heredia & Broderick Bownds
% dheredia@g.hmc.edu & brbownds@g.hmc.edu
% November 8, 2025

%% Question 1): Find A, mean value a0, estimate T0.

% Load data
load('HW9.mat'); % variable x

% Givens
Ts = 0.01;          % sec
Fs = 1 / Ts;        % Sampling Freq (Hz)
N = length (x);     % Number of samples
t = (0:N-1) * Ts;   % Time vector

% Plot HW9.mat signal
clc; clf;

figure(1);
plot(t, x, 'b', 'LineWidth', 2);
xlim([0,1]); ylim([-1,6])
grid on;
xlabel('Time (s)', 'FontSize', 20); ylabel('x(t)', 'FontSize', 20); title('Time-domain signal x(t)', 'FontSize', 20);

% Find a0
a0 = mean(x);

% Find A
LOW = median(x(x < a0));    % Approx LOW signal (stable held signal)
HIGH = median(x(x > a0));   % Approx HIGH signal (stable held signal)
A = HIGH - LOW;

% Note: Avoid using min/max since it assumes no noise spikes (unstable)

% Find T0 (Rising-Edge)
dx = diff(x);                           % Difference between consecutive samples
half_A = 0.5*A;                         % Half the pulse height
rising_edges = find(dx > half_A);       % Save index where rising edges occur

T_samples = diff(rising_edges);         % Difference between rising edges
T0_samples = median(T_samples);         % Take median to avoid outlier spikes
T0 = T0_samples * Ts;                   % Convert to sec

fprintf('Question 1): \na0 = %.3f, \nA = %.3f, \nT0 = %.3f\n', a0, A, T0);

% SOLUTION: a0 = 3.2381
% SOLUTION: A = 5.0000
% SOLUTION: T0 = 0.0900 s

%% Question 2): Plot magX vs f (Hz), plot magX vs k (-500 to 499).

% Defining parameters
x_samples = x(1:1000);              % First 1000 consecutive samples from x(t)
k = -500:499;                       % Given DFT index 
[X, f] = fdomain(x_samples, Fs);    % Call fdomain.m

% Plot magX vs f (Hz)
figure(2);
plot(f, abs(X), 'b', 'LineWidth', 2);
xlabel('Frequency (Hz)', 'FontSize', 20); ylabel('|X(f)|', 'FontSize', 20);
title('Magnitude Spectrum |X(f)| vs Frequency', 'FontSize', 20);
grid on;

% Plot magX vs k (-500 to 499)
figure(3);
plot(k, abs(X), 'r', 'LineWidth', 2);
xlabel('k', 'FontSize', 20); ylabel('|X(k)|', 'FontSize', 20);
title('Magnitude Spectrum |X(k)| vs Index', 'FontSize', 20);
grid on;

fprintf('\nQuestion 2): \nCheck plots in Git Repository\n');

% SOLUTION: Plot in Git Repository (Freq)
% SOLUTION: Plot in Git Repository (Index)

%% Question 3): Determine f0 from Question 2) and from f0 = 1/T0. Explain discrepancies.

% FFT of x_samples
X_samples = fft(x_samples);         % FFT of first 1000 consecutive samples          
magX = abs(X_samples);              % Magnitude of X (FFT)             

magX(1) = 0;                        % Omit DC value for fundamental frequency, not an avg value
maxVal = max(magX);                 % Find max magnitude and its index
idx = find(magX == maxVal, 1);      % Find index of first occurrence of max (signal peak)

% Convert index to freq
f0 = (idx - 1) * (Fs / length(x_samples));  
T0_fft = 1 / f0;

fprintf('\nQuestion 3): \nf0 = %.4f Hz, \nT0_fft = %.4f s, \nT0 = %.4f s\n', f0, T0_fft, T0);

% SOLUTION: f0 = 10.6000 Hz, 
% SOLUTION: T0_fft = 0.0943 s, 
% SOLUTION: T0 = 0.0900 s

% Explanation: The small, 0.0043 difference between T0_fft and T0 is due to
% the FFT providing a discrete approximation of the frequency of x(T). T0 
% is measured direclty from the time-domain rising edges, which reflects the
% actual pulse timing, whereas T0_fft is affected by the freq resolution
% of the FFT, sample length, and possible noise. T0 should be trusted more
% since it gives a more accurate measure of the period (T0 = 1 / f0).

%% Question 4): Determine spectrum of 1000pt Hanning, plot Plot magX vs f (Hz), plot magX vs k (-500 to 499).

% Defining paramters
hanning_win = hann(1000);           % Create Hanning window (1000pt)
x_hann = x_samples .* hanning_win;  % Define x_hann

% FFT of windowed signal
X_hann = fft(x_hann);              
magX_hann = abs(X_hann);            
X_shift = fftshift(X_hann);         % Shift zero-frequency to center

% Frequency vector
Fs = 1 / Ts;                        
N_hann = length(x_samples);              
f_hann = (0:N_hann-1)*(Fs/N_hann);  % Frequency axis (Hz)                      

% Plot magnitude vs frequency (Hz)
figure(4);
plot(f_hann, magX_hann, 'b', 'LineWidth', 2);
xlabel('Frequency (Hz)', 'FontSize', 20); ylabel('|X_{hann}(f)|', 'FontSize', 20);
title('Magnitude Spectrum |X_{hann}(f)| vs Frequency', 'FontSize', 20);
grid on;

% Plot magnitude vs index k (-500 to 499)
figure(5);

% Recall:
k = -500:499; 

plot(k, abs(X_shift), 'r', 'LineWidth', 2);
xlabel('k', 'FontSize', 20); ylabel('|X_{hann}(k)|', 'FontSize', 20);
title('Magnitude Spectrum |X_{hann}(k)| vs Index', 'FontSize', 20);
grid on;

fprintf('\nQuestion 4): \nCheck plots in Git Repository\n');

% SOLUTION: Plot in Git Repository (Freq)
% SOLUTION: Plot in Git Repository (Index)

%% Question 5): Calculate DTFS coefficients (a_k) for one period of x(t).

% Recall: A = 5.0000, a0 = 3.2381, T0 = 0.0900 s, Ts = 0.01 s

% Define parameters
d = a0 / A;                          % Duty cycle = a0 / A =  N1 / N0 
N1 = d * N0;                         % Number of samples per period when signal is ON
N0 = T0 / Ts;                        % Number of samples in one period

%%%%%%%%%%%%%%%%%%%%%
%                   %
%    UNFINISHED!    %
%                   %
%%%%%%%%%%%%%%%%%%%%% 
