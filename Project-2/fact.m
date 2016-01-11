function [x] = fact(n)

if n==1
    x = 1;
elseif n == 0
    x = 1;
else
    x = n * fact(n-1);
end
end