function y = fftx(x)
%FFTX Fast Finite Fourier Transform.
x = x(:);
n = length(x);
omega = exp(-2*pi*1i/n);

if rem(n,2) == 0

% Recursive divide and conquer.
k = (0:n/2-1)';
w = omega.^k;
u = fftx(x(1:2:n-1));
v = w.*fftx(x(2:2:n));
y = [u+v; u-v];
else
% The Fourier matrix.
j = 0:n-1;
k = j';
F = omega.^(k*j);
y = F*x;
end