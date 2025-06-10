%Load in Data 1) Using Filepath or 2) Drag File into Workspace
datapath = 'C:\Users\blah\Documents\blah\blah\test_scripts\my_data.mat';
%load(datapath)

% Each .mat file includes 
% theta: True angles of signal sources
% SNR: True SNRs of signal sources
% R_norm: Normalized Sample Covariance Matrix
% X_ml: AI/ML Surrogate Matrix Estimated from R
% N: Number of ULA Elements
% d: Element spacing of ULA as a multiple of wavelength
% az_extent: Extent in azimuth that network was trained under.


%Convert some data to double to avoid issues
az_extent = double(az_extent);
N = double(N);

%Create Array Manifold of Steering Vectors
az_angs = linspace(az_extent(1), az_extent(2), 500);
A=linear_dir_vec(N, d, az_angs);

%Calculate MUSIC Spectrum of R_norm
SssDim = min(length(theta), N-1);
[U,S,V]=svd(R_norm);
Un=U(:,SssDim+1:end);

for ii=1:length(az_angs)
    a=A(:,ii);
    Pmus(ii)=1/(a'*Un*Un'*a);
end

%Calculate Spatial Spectrum of X_ml
for ii=1:length(az_angs)
    a=A(:,ii);
    Psurr(ii)=1/(a'*X_ml*a);
end

% Plot both spectra
figure()
yyaxis left
hold on
plot(az_angs, 10*log10(abs(Pmus)))
for i = 1:length(theta)
    xline(theta(i), 'k--')
end
ylabel('MUSIC Pseudospectum (dB)')
yyaxis right
plot(az_angs, 10*log10(abs(Psurr)))
xlabel('Azimuth (deg)')
ylabel('Networks Spatial Spectrum (dB)')
grid on;grid minor


