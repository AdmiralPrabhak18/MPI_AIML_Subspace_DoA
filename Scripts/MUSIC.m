close all

N=6;               % Number of ULA Array Elements
d=0.5;             % Element Spacing (in terms of wavelength)
K=20;              % Number of Snapshots
az_angs=-40:.1:40; % Span of angles that are assessed 
theta=[3., -17, 0];% DoA of Sources, vector dim(1,M)
SNR=[40,40,40];    % Signal to Noise Ratio of Signals, vector dim(1,M)

SssDim = length(theta); % True Number of Sources, Used for Partitioning Eigenvector matrix
                        % Not available in practive. 

x = signal_gen (N, d, theta, SNR, K); % Array Data vector/matrix (depends on K)
R=x*x'/K;                             % Formation of Sample Covariance Matrix

[U,S,V]=svd(R);        % Eigendecomposition of R, U are eigenvectors, S are diag(eigenvalues)
Vs=U(:,1:SssDim);      % Signal Subspace
Vn=U(:,SssDim+1:end);  % Noise Subspace

R_inv = inv(R + eye(N)*8);
[~,S_inv, ~] = svd(R_inv);

eigValSpec=10*log10(diag(S)); % Eigenvalue Spectrum
eigValInvSpec=10*log10(diag(S_inv)); % Eigenvalue Spectrum

fprintf(1,'Signal Subspace Implemented= %2.0f, Number of Sources = %2.0f, Number of Elements =%2.0f\n',SssDim,length(theta), N);
fprintf(1,'Num.Samples = %2.0f, snr = %2.0f, %2.0f, %2.0f \n',K,SNR);
fprintf(1,'Source Angles = %2.1f, %2.1f, %2.1f \n',theta);
fprintf(1,'***** ***** *****\n');
fprintf(1,'Eigenvalues (dB) = [ %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f]\n',eigValSpec);

A=linear_dir_vec(N,d,az_angs); %Array Manifold of Steering Vectors
I=eye(N);
Pmus = [];
for ii=1:length(az_angs)
    a=A(:,ii); % Assess spectrum at the angle value corresponding to this particular steering vector
    Pmus(ii)=10*log10(abs(1/(a'*Vn*Vn'*a))); % Calculate MUSIC Spectrum Using Noise Subspace
    Pmus_signalSubSpace(ii)=10*log10(abs(1/(a'*(I-Vs*Vs')*a))); % Calculate Spectrum Using Signal Subspace
end
% Find Peaks in MUSIC Spectrum
[vals_mus, locs_mus] = findpeaks(Pmus);
peaks_mus = az_angs(locs_mus);
['Peaks From MUSIC: ',num2str(peaks_mus)]

%Plot EigenValue Spectrum
figure(1)
hold on
plot(eigValSpec, 'ko-','LineWidth',2)
plot(eigValInvSpec, 'go-','LineWidth',2)
title('Eigenvalue Spectrum')
xlabel('Eigenvalues')
ylabel('Magnitude (dB)')
grid on, zoom on, grid minor
ax = gca;
ax.XTick = unique(round(ax.XTick));


% Plot MUSIC Spectrum
figure(2)
hold on
plot(az_angs,(Pmus),'b', az_angs,(Pmus_signalSubSpace),'r--','LineWidth',2);
for i = 1:length(theta)
    xline(theta(i), '--')
end

grid on, zoom on
xlabel('Azimuth (deg)');
ylabel('PseudoSpectrum (dB)');
title('MUSIC Spectrum')
legend('Nss V_{n}','Sss V_{s}', 'Truth DoAs')



