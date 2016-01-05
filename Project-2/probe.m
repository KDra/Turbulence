flg = 1;
dt = 2e-3;
step = 1e-3;
counter = 1;
while flg
    dt = dt + step;
    flg = wave(dt);
    if ~flg && counter <= 4
        dt = dt - step;
        step = step/(10*counter);
        flg = 1;
        counter = counter + 1;
    end
end