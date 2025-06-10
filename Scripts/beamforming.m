% beamformULA
close all; clear all;
N = 6;                              % Number of ULA Array Elements
d = 0.5;                            % Element Spacing (in terms of wavelength)
theta=[-20, -10];                   % DoA of Sources, vector dim(1,M)
SNR=[30,30];                        % Signal to Noise Ratio of Signals, vector dim(1,M)
az_angs = -60:0.05:60;              % Span of angles that are assessed 
A = linear_dir_vec(N,d,az_angs);    % Array Manifold (set) of Steering Vector
x = signal_gen(N,d,theta,SNR,1);    % Single Snapshot of Array Data

% Estimation with beamformer
bf = x'*A;        % beamforming with data

% Plot Results
plot(az_angs,20*log10(abs(bf)),'linewidth',2);
for i=1:length(theta)
    xline(theta(i), '--','linewidth',1.5)
end
hold on;
legend('Beamformer', 'Truth DoA')
grid on; xlabel('Angle (degrees)'); ylabel('Magnitude(dB)');
title('Beamformer For Signal with ULA');
