flg = 1;        % Flag to stop the loop
dt = 1e-4;      % Initial step size
diff = 5e-3;    % Difference by which the step size is increased
counter = 1;    % Counter increases the resolution to more significant figures
while flg
    dt = dt + diff;
    % wave.m was modified to return a False value if the final energy
    % exceeds the initial value
    flg = wave(dt);     
    if ~flg && counter <= 4
        dt = dt - diff;             % Reverse previously taken step
        diff = diff/(10);           % Decrease difference
        flg = 1;                    % Set flag to True again
        counter = counter + 1;      % Increase the counter
    end
end