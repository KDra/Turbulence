N=10;
nu_range = NaN(2, N);
nu = linspace(1e-5, 0.01, N);
i=1;
cr = linspace(1e-5, 2, N);
for c = cr
    close all;
    flags = compNu(c, nu);
    pos1 = find(flags, 1, 'first');
    pos2 = find(flags, 1, 'last');
    if ~isempty(pos1) && ~isempty(pos2)
        nu_range(1, i) = nu(pos1);   % min
        nu_range(2, i) = nu(pos2);    % max
    end
    i=i+1;
end

