function [a, P, value] = pseudospettro(w, N, G)
a = exp(-j*[0:N-1]'*w);     %steering vector relativo ad w
P = 1./(a'*(G)*(G)'*a);     %pseudospettro relativo a a(w)
value = abs(P);
