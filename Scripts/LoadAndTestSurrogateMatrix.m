%Load in Data 1) Using Filepath or 2) Drag File into Workspace
datapath = 'C:\Users\blah\Documents\blah\blah\test_scripts\my_data.mat';
%load(datapath)

% Each .mat file includes 
% theta: True angles of signal sources
% SNR: True SNRs of signal sources
% K: Number of Snapshots
% x: Array Data vector/matrix (dim NxK)
% R_norm: Normalized Sample Covariance Matrix
% X_ml: AI/ML Surrogate Matrix Estimated from R
% N: Number of ULA Elements
% d: Element spacing of ULA as a multiple of wavelength
% az_extent: Extent in azimuth that network was trained under.


%Convert some data to double to avoid issues
az_extent = double(az_extent);
N = double(N);
K = double(K);

%Create Array Manifold of Steering Vectors
az_angs = linspace(az_extent(1), az_extent(2), 500);
A=linear_dir_vec(N, d, az_angs);

%Calculate MUSIC Spectrum of R_norm
SssDim = min(length(theta), N-1); %Dimension of Signal Subspace
[U,S,V]=svd(R_norm);  % Eigen Decomposition
Un=U(:,SssDim+1:end); % Form Noise Subspace
Pmus = [];
for ii=1:length(az_angs)
    a=A(:,ii);  %Select a steering vector
    Pmus(ii)=10*log10(abs(1/(a'*Un*Un'*a)));  % Calculate MUSIC spectra
end
% Find Peaks in MUSIC Spectrum
[vals_mus, locs_mus] = findpeaks(Pmus);
peaks_mus = az_angs(locs_mus);

%Calculate Spatial Spectrum of X_ml
Pml = [];
for ii=1:length(az_angs)
    a=A(:,ii);
    Pml(ii)=10*log10(abs(1/(a'*X_ml*a))); % Calculate spatial spectra of Xml
end
% Find peaks estimated from Xml
[vals_ml, locs_ml] = findpeaks(Pml);
peaks_ml = az_angs(locs_ml);

% Plot both spectra
figure()
yyaxis left

hold on
scatter(peaks_mus, vals_mus + 5, 'vb')
plot(az_angs, Pmus, 'b')
for i = 1:length(theta)
    xline(theta(i), 'k--')
end
ylabel('MUSIC Pseudospectum (dB)')
yyaxis right
plot(az_angs, Pml, 'r')
scatter(peaks_ml, vals_ml + 5, 'vr')
xlabel('Azimuth (deg)')
ylabel('Networks Spatial Spectrum (dB)')
grid on;grid minor

disp(['Peaks From MUSIC: ',num2str(peaks_mus)])
disp(['Peaks From Neural Networks: ',num2str(peaks_ml)])


