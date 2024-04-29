day      = 60*60*24; % Day length (s).
tmax     = day * 10; % Duration of the simulation (s).
clockmax = 400 ;% Number of time steps.
dt = tmax/clockmax ;% Calculates the duration of each time step.


%% 2.5 births per death --> 1/2.5 deaths per birth 
a           = 100/day;             % infectivity 
b           = 0.1/day;             % recovery rate 
betaHealth  = .001/day;     % birthrate for healthy
betaInf     = betaHealth* (1/4)  ;% birthrate for ill 
deltaHealth = 1/2.5 * betaHealth; 
deltaInf    = deltaHealth* (5);      % death rate for infected individuals 
rv          = .00/day;

%% Note that betaInf < betaHealth | deltaInf > deltaHealth | 
% ... betaHealth > deltaHealth | deltaInf > betaInf

N = 1000    ;% Total population
I = 100     ;% Infected
S = N - I   ;% Susceptible 
R = 0       ;% Recovered
D = 0       ;% Deceased
V = 0       ;% Vaccinated

tsave = zeros(1,clockmax);
Ssave = zeros(1,clockmax);
Isave = zeros(1,clockmax);
Rsave = zeros(1,clockmax);
Dsave = zeros(1,clockmax);
Vsave = zeros(1,clockmax);
Nsave = zeros(1,clockmax);


figure;
hold on;

hS = plot(0, S, 'g', 'linewidth', 2);
hI = plot(0, I, 'r', 'linewidth', 2);
hR = plot(0, R, 'b', 'linewidth', 2);
hD = plot(0, D, 'k', 'linewidth', 2);
hV = plot(0, V, 'y', 'linewidth', 2);
hN = plot(0, N, 'm', 'linewidth', 2);


drawnow

legend({'S','I','R', 'D', 'V'},'Location','northeast')

axis([0, tmax, 0, 1.02])
%% ds/dt = -a(ptrans)*S + betaHealth*S + betaInf*I - deltaHealth*S
%% di/dt = a(ptrans)*S - betaInf*I - deltaInf*S
%% N = S + I + R 
for clock=1:clockmax
    t = clock*dt; % Updates current time

    ptrans = I/N;
    % Calculates new cases of I, R and D
    if S > 0
        newI    = dt * a * ptrans;
        Sbirths = dt * betaHealth * S;
        Sdie    = dt * deltaHealth* S;
        Ibirths = dt * betaInf*I;
        newV    = dt * rv * S;
    else
        newI = 0;
        newV = 0;
    end
    if I > 0
        newR = dt*b*I;
        newD = dt*deltaInf*I;
    else
        newR = 0;
        newD = 0;
    end
    
    % Calculate final values of variables
    S = S - newI + Sbirths + Ibirths - Sdie - newV;
    I = I + newI - newR    - newD;
    R = R + newR;
    D = D + newD;
    V = V + newV;
    N = S + I + R + V;
    
    % Update tsave, Ssave, Isave, Rsave, Dsave
    tsave(clock) = t; 
    Ssave(clock) = S/N;
    Isave(clock) = I/N;
    Rsave(clock) = R/N;
    Dsave(clock) = D/N;
    Vsave(clock) = V/N;
    % Nsave(clock) = N;

    check = (S+I+R+V) / N
    % Update the plots
    set(hS, 'XData', tsave(1:clock), 'YData', Ssave(1:clock));
    set(hI, 'XData', tsave(1:clock), 'YData', Isave(1:clock));
    set(hR, 'XData', tsave(1:clock), 'YData', Rsave(1:clock));
    set(hD, 'XData', tsave(1:clock), 'YData', Dsave(1:clock));
    set(hV, 'XData', tsave(1:clock), 'YData', Vsave(1:clock));

    % set(hN, 'XData', tsave(1:clock), 'YData', Nsave(1:clock));

    
    drawnow update; % Update the plot
end