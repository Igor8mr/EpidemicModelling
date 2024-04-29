day      = 60*60*24; % Day length (s).
tmax     = day * 100; % Duration of the simulation (s).
clockmax = 400 ;% Number of time steps.
dt = tmax/clockmax ;% Calculates the duration of each time step.

%% Model Parameters
A           = 1/day;  % infectivity 
B           = 0.01/day;  % recovery rate 

a           = [A, A/2, 0];
b           = [B, B/2, B];
ra          = 0.5;         % reinfection multiplier

betaH       = 0.001/day;   % birthrate for healthy
betaI       = betaH * (1/4); % birthrate for ill 

deltaH      = betaH;       % Death rate for healthy individuals
deltaI      = [deltaH * 5, deltaH * 5 /2,  deltaH * 5]; % Death rate for infected individuals

vr          = 1/day;       % Recovery rate
qr          = 1/day;       % Quarantine rate

%% Initial Conditions
N =     [1000,        0,     0] ; % Total population
I =     [100,         0,     0] ; % Infected
S =     [N(1)-I(1),   0,     0] ; % Susceptible 
R =     [0,           0,     0] ; % Recovered
D =     [0,           0,     0] ; % Total Deceased

%% Initialization for Plotting
tsave = zeros(clockmax, 1);
Ssave = zeros(clockmax, 3);
Isave = zeros(clockmax, 3);
Rsave = zeros(clockmax, 3);
Dsave = zeros(clockmax, 3);
Nsave = zeros(clockmax, 3);

%% Create the figure and subplots
figure;

subplot(2,2,1);

hold on

hS1 = plot(tsave(1:clockmax), sum(Ssave(1:clockmax)));
hI1 = plot(tsave(1:clockmax), sum(Isave(1:clockmax)));
hR1 = plot(tsave(1:clockmax), sum(Rsave(1:clockmax)));
hD1 = plot(tsave(1:clockmax), sum(Dsave(1:clockmax)));

legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax, 0, 1.02])
title('Subplot 1')

subplot(2,2,2);

hold on

hS2 = plot(tsave(1:clockmax), Ssave(1:clockmax, 1));
hI2 = plot(tsave(1:clockmax), Isave(1:clockmax, 1));
hR2 = plot(tsave(1:clockmax), Rsave(1:clockmax, 1));
hD2 = plot(tsave(1:clockmax), Dsave(1:clockmax, 1));

legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax, 0, 1.02])


subplot(2,2,3);

hold on

hS3 = plot(tsave(1:clockmax), Ssave(1:clockmax, 2));
hI3 = plot(tsave(1:clockmax), Isave(1:clockmax, 2));
hR3 = plot(tsave(1:clockmax), Rsave(1:clockmax, 2));
hD3 = plot(tsave(1:clockmax), Dsave(1:clockmax, 2));

legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax, 0, 1.02])


subplot(2,2,4);

hold on

hS4 = plot(tsave(1:clockmax), Ssave(1:clockmax, 3));
hI4 = plot(tsave(1:clockmax), Isave(1:clockmax, 3));
hR4 = plot(tsave(1:clockmax), Rsave(1:clockmax, 3));
hD4 = plot(tsave(1:clockmax), Dsave(1:clockmax, 3));

legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax, 0, 1.02])
title('Subplot 1')

drawnow update;

%% Main Simulation Loop
for clock = 1:clockmax
    t = clock*dt; % Updates current time

    % Calculating populational changes
    ptrans = (I(1) + I(2)) / (N(1) + N(2));

    Sbirths = dt * (betaH * (sum(S)+sum(R)) + betaI * sum(I));
    Sinf = dt * ptrans * a .* S;
    Sdie = dt * deltaH * S;

    Idie = dt * deltaI .* I;

    Rnew = dt * b .* I;
    Rinf = dt * ptrans * ra * a .* R;
    Rdie = dt * deltaH * R;
    
    % Calculating final values of variables
    S = S + [Sbirths, 0, 0] - Sdie - Sinf;
    I = I + Sinf + Rinf - Idie - Rnew;
    R = R + Rnew - Rinf - Rdie;
    D = D + Sdie + Idie + Rdie;

    S(1) = S(1) + S(1) * (qr+vr) * dt;
    S(3) = S(3) + S(1) * qr * dt - S(3) * vr * dt;
    S(2) = S(2) + S(1) * vr * dt + S(3) * vr * dt;
    
    I(1) = I(1) - I(1) * qr * dt;
    I(3) = I(3) + I(1) * qr * dt;

    R(1) = R(1) + R(1) * (qr+vr) * dt;
    R(3) = R(3) + R(1) * (ptrans-0.1) * qr * dt - R(3) * vr * dt;
    R(2) = R(2) + R(1) * vr * dt + R(3) * vr * dt;

    N = S + R + I;
    
    % Update tsave, Ssave, Isave, Rsave, Dsave
    tsave(clock) = t; 
    Ssave(clock, :) = S ./ N;
    Isave(clock, :) = I ./ N;
    Rsave(clock, :) = R ./ N;
    Dsave(clock, :) = D ./ N;

    % Update the plots in the first subplot
    subplot(2,2,1)
    set(hS1, 'XData', tsave(1:clock), 'YData', sum(Ssave(1:clock)));
    set(hI1, 'XData', tsave(1:clock), 'YData', sum(Isave(1:clock)));
    set(hR1, 'XData', tsave(1:clock), 'YData', sum(Rsave(1:clock)));
    set(hD1, 'XData', tsave(1:clock), 'YData', sum(Dsave(1:clock)));

    % Update the plots in the second subplot
    subplot(2,2,2)
    set(hS2, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 1));
    set(hI2, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 1));
    set(hR2, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 1));
    set(hD2, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 1));

     % Update the plots in the first subplot
    subplot(2,2,3)
    set(hS1, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 2));
    set(hI1, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 2));
    set(hR1, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 2));
    set(hD1, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 2));

    % Update the plots in the second subplot
    subplot(2,2,4)
    set(hS2, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 3));
    set(hI2, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 3));
    set(hR2, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 3));
    set(hD2, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 3));
    
    drawnow; % Update the plot
end
