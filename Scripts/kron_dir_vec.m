function y_kron = kron_dir_vec (N, d, theta)
% assumes electronic steering of zero degrees (boresight)
% INPUTS
% N = number of elements
% d = element spacing (wavelengths)
% theta = angle (degrees) (this can be a vector)
% plotFlag = boolean to determine whether to plot pattern

% OUTPUTS
% y = radiation pattern vector(ne x len(theta)) (complex voltage)
    y = linear_dir_vec(N, d, theta);
    N = size(y, 1);
    num_ang = size(y, 2);

    y_kron = zeros(N.^2, num_ang);
    y_kron = complex(y_kron);
    for i = 1:num_ang
        a = y(:, i);
        y_kron(:,i) = kron(a, conj(a));
    end
end


