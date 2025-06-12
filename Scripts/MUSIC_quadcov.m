close all

N=4;               % Number of ULA Array Elements
d=0.5;             % Element Spacing (in terms of wavelength)
K=10000;             % Number of Snapshots
az_angs=-90:.1:90; % Span of angles that are assessed 
theta=[-25, -10, 10, 25];% DoA of Sources, vector dim(1,M)
SNR=[40, 40, 40, 40];    % Signal to Noise Ratio of Signals, vector dim(1,M) 

[x, Qc, QcFB] = signal_gen_quadcov (N, d, theta, SNR, K); % Array Data vector/matrix (depends on K)

% Obtain Noise Subspace Projection of R
SssDim_R = min(length(theta), N-1); % True Number of Sources, Used for Partitioning Eigenvector matrix
                        % Not available in practive.
R=x*x'/K;              % Formation of Sample Covariance Matrix
[U_R,S_R,V_R]=svd(R);        % Eigendecomposition of R, U are eigenvectors, S are diag(eigenvalues)
Vs_R=U_R(:,1:SssDim_R);      % Signal Subspace of R
Vn_R=U_R(:,SssDim_R+1:end);  % Noise Subspace of R

eigValSpec_R=10*log10(diag(S_R)); % Eigenvalue Spectrum

% Calculate MUSIC Spatial Spectra
A_cov=linear_dir_vec(N,d,az_angs); %Array Manifold of Steering Vectors
Pmus = [];
for ii=1:length(az_angs)
    a=A_cov(:,ii); % Assess spectrum at the angle value corresponding to this particular steering vector
    Pmus(ii)=10*log10(abs(1/(a'*Vn_R*Vn_R'*a))); % Calculate MUSIC Spectrum Using Noise Subspace
end
% Find Peaks in MUSIC Spectrum
[vals_mus, locs_mus] = findpeaks(Pmus);
peaks_mus = az_angs(locs_mus);
['Peaks From MUSIC on Sample Covariance: ',num2str(peaks_mus)]

% Obtain Noise Subspace Projection of Qc
SssDim_quad = length(theta); % True Number of Sources, Used for Partitioning Eigenvector matrix
                           % Not available in practive.
[U_quad,S_quad,V_quad]=svd(Qc);        % Eigendecomposition of Qc, U are eigenvectors, S are diag(eigenvalues)
Vs_quad=U_quad(:,1:SssDim_quad);      % Signal Subspace of Qc
Vn_quad=U_quad(:,SssDim_quad+1:end);  % Noise Subspace of Qc

eigValSpec_quad=10*log10(diag(S_quad)); % Eigenvalue Spectrum

% Calculate MUSIC Spatial Spectra for QUadricovariance Matrix
A_quadcov=kron_dir_vec(N,d,az_angs); %Array Manifold of Steering Vectors
Pquad = [];
for ii=1:length(az_angs)
    a=A_quadcov(:,ii); % Assess spectrum at the angle value corresponding to this particular steering vector
    Pquad(ii)=10*log10(abs(1/(a'*Vn_quad*Vn_quad'*a))); % Calculate MUSIC Spectrum Using Noise Subspace
end
% Find Peaks in MUSIC Spectrum
[vals_quad, locs_quad] = findpeaks(Pquad);
peaks_quad = az_angs(locs_quad);
['Peaks From MUSIC on Quadricovariance : ',num2str(peaks_quad)]

%Plot EigenValue Spectra
figure(1)
hold on
plot(eigValSpec_R, 'bo-','LineWidth',2)
plot(eigValSpec_quad, 'ro-','LineWidth',2)
title('Eigenvalue Spectrum')
xlabel('Eigenvalues')
ylabel('Magnitude (dB)')
grid on, zoom on, grid minor
ax = gca;
ax.XTick = unique(round(ax.XTick));
legend('Correlation','Quadricovariance')


% Plot Spatial Spectra
figure(2)
hold on
plot(az_angs,(Pmus),'b','LineWidth',2);
plot(az_angs,(Pquad),'r','LineWidth',1);
for i = 1:length(theta)
    xline(theta(i), '--')
end
grid on, zoom on
xlabel('Azimuth (deg)');
ylabel('PseudoSpectrum (dB)');
title('MUSIC Spectrum')
legend('MUSIC (Correlation)','MUSIC (Quadricovariance)', 'Truth DoAs')


