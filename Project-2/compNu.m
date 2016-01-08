function [flags] = compNu(c, nu)
N = size(nu, 2);
flags = ones(1,N);
for i = 1:N
    flags(i) = wave(c, nu(i));
end
end