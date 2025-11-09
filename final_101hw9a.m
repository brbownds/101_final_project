% E101: HW 9a
% Sebastian Heredia & Broderick Bownds
% dheredia@g.hmc.edu & brbownds@g.hmc.edu
% November 8, 2025

% Question 1): Find A, mean value a0, estimate T0

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
plot(t, x, 'b');
xlim([0,1]); ylim([-1,6])
grid on;
xlabel('Time (s)'); ylabel('x(t)'); title('Time-domain signal x(t)');

% Find a0
a0 = mean(x);

% Find A
LOW = median(x(x < a0));    % Approx LOW signal (stable held signal)
HIGH = median(x(x > a0));   % Approx HIGH signal (stable held signal)
                            % Using min/max assumes no noise spikes (unstable)

A = HIGH - LOW;

% Find T0 (Rising-Edge)
dx = diff(x);                           % Difference between consecutive samples
half_A = 0.5*A;                         % Half the pulse height
rising_edges = find(dx > half_A);       % Save index where rising edges occur

T_samples = diff(rising_edges);         % Difference between rising edges
T0_samples = median(T_samples);         % Take median to avoid outlier spikes
T0 = T0_samples * Ts;                   % Convert to sec

fprintf('Question 1) : a0 = %.3f, A = %.3f, T0 = %.3f\n', a0, A, T0);
% SOLUTION: a0 = 3.2381
% SOLUTION: A = 5.0000
% SOLUTION: T0 = 0.0900 s

% ///////////////////////////////////////////////////////////////////////////////////// %

% Question 2): Choose 1000 data points from "HW9.mat" to perform the
% FFT using fdomain.m to plot magnitude vs frequency (Hz). Also 

x_section = x(1:1000); % data points

[X, f] = fdomain(x_section, Fs); % perform the FFT where X gives the complex 
                                % FFT result and f is the frequency vector
% Plot Magnitude vs Frequency
figure (2);
plot(f, abs(X), 'r');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
title('Magnitude Specturm vs. Frequency');

% Plot Magnitude vs. DFT index where k = [-500, 499]
X_shifted = fftshift(X);
k = -500:499;

figure(3);
plot(k, abs(X_shifted));
xlabel('DFT index (k)');
ylabel('|X(k)|');
title('Magnitude Spectrum vs. DFT index where k = [-500, 499]');
