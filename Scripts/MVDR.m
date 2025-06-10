close all

N=6;               % Number of ULA Array Elements
d=0.5;             % Element Spacing (in terms of wavelength)
K=20;              % Number of Snapshots
az_angs=-40:.1:40; % Span of angles that are assessed 
theta=[3., -17, 0];% DoA of Sources, vector dim(1,M)
SNR=[40,40,40];    % Signal to Noise Ratio of Signals, vector dim(1,M)


x = signal_gen (N, d, theta, SNR, K); % Array Data vector/matrix (depends on K)
R=x*x'/K;                             % Formation of Sample Covariance Matrix

A=linear_dir_vec(N,d,az_angs); %Array Manifold of Steering Vectors
I=eye(N);
for ii=1:length(az_angs)
    a=A(:,ii); % Assess spectrum at the angle value corresponding to this particular steering vector
    CAPON(ii)=1/(a'*inv(R)*a); % Perform inner product of steering vector against inverse of R
end

% Plot MVDR Spectrum
figure(1)
hold on
plot(az_angs, 10*log10(CAPON) ,'LineWidth',2)
for i = 1:length(theta)
    xline(theta(i), '--')
end

grid on, zoom on
xlabel('Azimuth (deg)');
ylabel('PseudoSpectrum (dB)');
title('MVDR Spectrum')
legend( 'MVDR', 'Truth DoAs')


