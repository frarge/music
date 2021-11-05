Ny1 = 91:20:291;    %lunghezza segnale osservato
N = 11;     %lunghezza segnale da stimare
K = 4;      %numero componenti tonali
w = [-0.3 -0.25 0.125 0.325]*2*pi;
alpha = [2 5 3 4];
ALPHA = [];
phi = rand(1,K)*2*pi; 
sigma_noise = 1;
NFFT = 1024;

A_hat = [];
F_hat = [];

for Ny = Ny1
    y = 0; n = 0:Ny-1;
    for k = 1:K
        y = y+alpha(k)*exp(j*(w(k)*n+phi(k)));
    end
    %generazione del segnale
    y = y + sigma_noise*randn(1,Ny);
    [P, F, w_hat] = music(y,N,K,NFFT);
    window = [ones(1,Ny)];
    pyy = periodogram(y,window,2*pi*F)';
    figure, plot(F, log(abs(P)),F,log(abs(pyy)));   %scala logaritmica
    xlabel('$F$', 'interpreter', 'latex')
    ylabel('$|P|$', 'interpreter', 'latex')
    legend('MUSIC','PERIODOGRAM');
    title({'Stima frequenze con Ny=',Ny})
    %stima ampiezze con minimi quadrati
    B = [];
    for l = 1:K
        b = exp(j*w_hat(l)).^(0:Ny-1);
        B = [B b.'];
    end
    alpha_hat = pinv(B)*y(:);
    A_hat = [A_hat alpha_hat];
    disp('stima ampiezze:')
    disp(abs(alpha_hat'))
    disp('stima frequenze:')
    disp(w_hat/(2*pi))
    F_hat = [F_hat (w_hat/(2*pi))'];
    ALPHA = [ALPHA alpha'];
end
%plot risultati
figure, plot(Ny1, abs(A_hat),Ny1,ALPHA)
xlabel('$Ny$', 'interpreter', 'latex')
ylabel('$\hat{\alpha}$', 'interpreter', 'latex')
legend('\alpha_1','\alpha_2','\alpha_3','\alpha_4')
title(['Stima Ampiezze'])
figure, plot(Ny1, F_hat)
xlabel('$Ny$', 'interpreter', 'latex')
ylabel('$\hat{w}$', 'interpreter', 'latex')
legend('f_1','f_2','f_3','f_4')
title(['Stima Frequenze'])
