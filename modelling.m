day = 60*60*24 % Day length (s).
tmax = day * 10 % Duration of the simulation (s).
clockmax = 400 % Number of time steps.
dt = tmax/clockmax % Calculates the duration of each time step.

a = 500/day
b = 0.5/day
c = 0.1/day

N = 1000    % Total population
I = 100     % Infected
S = N - I   % Susceptible 
R = 0       % Recovered
D = 0       % Deceased

tsave = zeros(1,clockmax);
Ssave = zeros(1,clockmax);
Isave = zeros(1,clockmax);
Rsave = zeros(1,clockmax);
Dsave = zeros(1,clockmax);

figure;
hold on;

hS = plot(0, S, 'g', 'linewidth', 2);
hI = plot(0, I, 'r', 'linewidth', 2);
hR = plot(0, R, 'b', 'linewidth', 2);
hD = plot(0, D, 'k', 'linewidth', 2);

legend({'S','I','R', 'D'},'Location','northeast')

axis([0, tmax, 0, 1.05 * N])

for clock=1:clockmax
    t = clock*dt; % Updates current time

    ptrans = I/N;

    % Calculates new cases of I, R and D
    if S > 0
        newI = dt*a*ptrans;
    else
        newI = 0;
    end
    if I > 0
        newR = dt*b*I;
        newD = dt*c*I;
    else
        newR = 0;
        newD = 0;
    end
    
    % Calculate final values of variables
    S = S - newI
    I = I + newI - newR - newD
    R = R + newR
    D = D + newD
    
    % Update tsave, Ssave, Isave, Rsave, Dsave
    tsave(clock) = t; 
    Ssave(clock) = S;
    Isave(clock) = I;
    Rsave(clock) = R;
    Dsave(clock) = D;
    
    % Update the plots
    set(hS, 'XData', tsave(1:clock), 'YData', Ssave(1:clock));
    set(hI, 'XData', tsave(1:clock), 'YData', Isave(1:clock));
    set(hR, 'XData', tsave(1:clock), 'YData', Rsave(1:clock));
    set(hD, 'XData', tsave(1:clock), 'YData', Dsave(1:clock));
    
    drawnow update; % Update the plot
end