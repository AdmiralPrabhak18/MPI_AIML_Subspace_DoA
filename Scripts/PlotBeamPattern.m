close all
N_1 = 4;    % Number of Antenna Elements in ULA 1
N_2 = 12;   % Number of Antenna Elements in ULA 2
d = 0.5;    % Element Spacing (in terms of wavelength)
az_angs = linspace(-180,180, 1000); % Span of angles that are assessed
eps=1e-6;   % Remove log(0)

A_1 = linear_dir_vec(N_1, d, az_angs); %Array Manifold of Steering Vectors
boresight_1 = linear_dir_vec(N_1, d, 0); % Steering vector at boresight (az=0)
                                         % Equivalent to forming sum beam
                                         % (summation of all elements)
AntPattern_1 = abs(boresight_1' * A_1 / sqrt(N_1)); % Form antenna response across angle extent

% Repeat for Array 2
A_2 = linear_dir_vec(N_2, d, az_angs); 
boresight_2 = linear_dir_vec(N_2, d, 0);
AntPattern_2 = abs(boresight_2' * A_2 / sqrt(N_2));

% Plot Antenna Patterns
figure()
hold on
plot(az_angs, 20*log10(AntPattern_1 / max(AntPattern_1)),'LineWidth',2)
plot(az_angs, 20*log10(AntPattern_2/ max(AntPattern_2)),'LineWidth',2)
ylim([-60,0]);
grid on; grid minor 
xlabel('Azimuth (deg)')
ylabel('Normalized Response (dB)')
legend({'4 Element ULA', '12 Element ULA'})

% Plot Antenna Patterns (Polar Plots)
figure()
polarplot(az_angs*(2*pi) / 360,  AntPattern_1 ,'LineWidth',2)
hold on 
polarplot(az_angs*(2*pi) / 360,  AntPattern_2 / min(AntPattern_2),'LineWidth',2)
hold off
legend({'4 Element ULA', '12 Element ULA'})
grid on; grid minor 
%%