N=5;
nu_range = NaN(2, N);
nu = linspace(0.0066, 0.007, N);
i=1;
for c = linspace(1.7, 1.8, N)
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

