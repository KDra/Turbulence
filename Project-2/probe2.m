
function [lamda, mu] = probe2()
N=30;      % No of iterations
nu_range = NaN(2, N);
nu = linspace(1e-5, 0.01, N);
i=1;
cr = linspace(1e-5, 4, N);
% Compute for each value of c
for c = cr
    flags = compNu(c, nu);
    pos1 = find(flags, 1, 'first');
    pos2 = find(flags, 1, 'last');
    if ~isempty(pos1) && ~isempty(pos2)
        nu_range(1, i) = nu(pos1);   % min
        nu_range(2, i) = nu(pos2);    % max
    end
    i=i+1;
end
% Obtain CFL values
dt = 0.003;
dx = 0.0073;
lamda = nu_range * dt/dx^2;
mu = cr * dt/dx;
plot(cr, lamda(1, :), cr, lamda(2, :))
end

function [flags] = compNu(c, nu)
% See if the function is stable for the given c and nu values
% c     float       a wavespeed value
% nu    vector      spaced viscosity values
% flags vector      1 - stable, 0 - unstable
N = size(nu, 2);
flags = ones(1,N);
for i = 1:N
    flags(i) = wave(c, nu(i));
end
end