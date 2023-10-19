function [extinct,time] = stochastic_simulation(x_birth,x_death, ...
    y_birth,y_death,volume,x,y)
%STOCHASTIC_SIMULATION Simulates a trajectory of the continuous time random
%walk with Lotka-Volterra transition rates. Outputs TRUE if y = 0 and the
%total simulation time.

time = 0;
extinct = false;
marginal = false;
minimum = false;

while ~extinct && ~marginal

    % Compute transition rates
    r1 = x_birth * x * y / volume;
    r2 = x_death * x;
    r3 = y_birth * y;
    r4 = y_death * x * y / volume;
    r0 = r1 + r2 + r3 + r4;

    % Update reaction time
    time = time + log(1 / rand) / r0;

    % Select transition and update
    r = rand;
    if r < r1 / r0
        x = x + 1;
    elseif r < (r1 + r2) / r0
        x = x - 1;
    elseif r < (r1 + r2 + r3) / r0
        y = y + 1;
    else
        y = y - 1;
    end

    % Determine if minimum reached
    if x < y_birth * volume / y_death && y < x_death * volume / x_birth
        minimum = true;
    end

    % Determine if in marginal state
    if (minimum && y > x_death * volume / x_birth) || x == 0
        marginal = true;
    end

    % Check for extinction event
    if y == 0
        extinct = true;
    end

end

