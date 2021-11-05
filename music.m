function [P, F, w_hat] = MUSIC(y,N,K,NFFT)
R = covar(y,N);
[V, lambda] = eig(R);   %calcolo autovettori di Ryy
[s, ind] = sort(diag(lambda));
G = V(:,ind((1:N-K)));      %estraggo autovettori appartenenti al sottospazio di rumore
P = [];
F = -0.5:1/NFFT:0.5;    %campionamento fitto fra -pi e pi
w_hat = [];     %contiene stima pulsazioni

%calcolo dei punti di massimo e pseudo-spettro
while(size(w_hat,2) < K)
    for w = 2*pi*F
        [a, L, value] = pseudospettro(w, N, G);
        P = [P L];
        w_prec = w-((1/NFFT)*2*pi);
        w_next = w+((1/NFFT)*2*pi);
        [a_prec,L_prec,v_prec] = pseudospettro(w_prec,N,G);
        [a_next,L_next,v_next] = pseudospettro(w_next,N,G);
        diff1 = value-v_prec;
        diff2 = v_next-value;
        if(diff1*diff2 < 0 && value>v_prec && diff1>0.0001 && size(w_hat,2) < K)
            w_hat = [w_hat w];
        end
    end
    if(size(w_hat) < K)
        P = [];
    end
end
return