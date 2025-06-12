function [x, Qc, QcFB] = signal_gen_quadcov(N, d, theta, SNR, K)
x = signal_gen(N, d, theta, SNR, K, true);

shape = zeros(N,N,N,N);
C = complex(shape,0);
shape = zeros(N.^2,N.^2);
Qc = complex(shape,0);

for i = 1:N
    for j = 1:N
        for k = 1:N
            for l = 1:N
                index1 = N*(i-1) + k;
                index2 = N*(l-1) + j;
                
                c1 = (sum(x(i,:).*conj(x(j,:)).*conj(x(k,:)).*x(l,:))) / K;
                c2 = ( sum(x(i,:).*conj(x(j,:))) .* sum(conj(x(k,:)).*x(l,:)) ) / K.^2;
                c3 = ( sum(x(i,:).*conj(x(k,:))) .* sum(conj(x(j,:)).*x(l,:)) ) / K.^2;
                c4 = ( sum(x(i,:).*x(l,:)) .* sum(conj(x(j,:)).*conj(x(k,:))) ) / K.^2;
                
                C(i,j,k,l) = c1-c2-c3-c4;
                
                Qc(index1, index2) = C(i,j,k,l);
                
            end
        end
    end
end
            
anti = fliplr(diag(ones(1, N.^2)));
QcB = anti * conj(Qc) * anti;
QcFB = (Qc + QcB)/2;
%%
end