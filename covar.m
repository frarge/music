function R = covar(y,N)
y = y(:);
Ny = numel(y);
Y = [];
for k=N:Ny
    yn = y(k:-1:k-N+1);
    Y = [Y yn];
end
R = (Y*Y')/Ny;
return