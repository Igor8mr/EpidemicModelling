day      = 60*60*24; % Day length (s).
tmax     = day * 365; % Duration of the simulation (s).
clockmax = 400 ;% Number of time steps.
dt = tmax/clockmax ;% Calculates the duration of each time step.


%% 2.5 births per death --> 1/2.5 deaths per birth 
A           = 1/day;  %  infectivity 
B           = 1/day;  %  recovery rate 

a           = [A, A/2, 0];
b           = [B, B/2, B];
ra          = 0.5;         % reinfection multiplier

betaH  =      0.001/day;         % birthrate for healthy
betaI     =   betaH * (1/4);     % birthrate for ill 

deltaH      = betaH;    % Death rate for healthy individuals
deltaI      = [deltaH * 5, deltaH * 5 /2,  deltaH * 5]; % Death rate for infected individuals

vr          = 1/day;

%% Note that betaI < betaH | deltaI > deltaH | 
% ... betaH > deltaH | deltaI > betaI

%       Non-vaccinated
%                     Vaccinated
%                            Quarentined and non-vaccinated
N =     [1000,        0,     0] ; % Total population
I =     [100,         0,     0] ; % Infected
S =     [N(1)-I(1),   0,     0] ; % Susceptible 
R =     [0,           0,     0] ; % Recovered
D =     [0,           0,     0] ; % Total Deceased

tsave = zeros(1,clockmax);
Ssave = zeros(1,clockmax);
Isave = zeros(1,clockmax);
Rsave = zeros(1,clockmax);
Dsave = zeros(1,clockmax);
Vsave = zeros(1,clockmax);
Nsave = zeros(1,clockmax);

figure;
hold on;

hS = plot(0, sum(S), 'g', 'linewidth', 2);
hI = plot(0, sum(I), 'r', 'linewidth', 2);
hR = plot(0, sum(D), 'b', 'linewidth', 2);
hD = plot(0, sum(D), 'k', 'linewidth', 2);
hV = plot(0, sum(V), 'y', 'linewidth', 2);
hN = plot(0, sum(N), 'm', 'linewidth', 2);

drawnow

legend({'S','I','R', 'D', 'V'},'Location','northeast')

axis([0, tmax, 0, 1.02])
%% ds/dt = -a(ptrans)*S + betaH*S + betaI*I - deltaH*S
%% di/dt = a(ptrans)*S - betaI*I - deltaI*S
%% N = S + I + R 
for clock=1:clockmax
    t = clock*dt; % Updates current time

    %% Calculating populational changes
    ptrans =    (I(1) + I(2)) / (N(1) + N(2))

    Sbirths =   dt * (betaH * (sum(S)+sum(R)) + betaI * sum(I));
    Sinf =      dt * ptrans * a .* S;
    Sdie =      dt * deltaH * S;

    Idie =      dt * deltaI .* I;

    Rnew =      dt * b .* I;
    Rinf =      dt * ptrans * ra * a .* R;
    Rdie =      dt * deltaH * R;
    
    %% Calculating final values of variables
    S = S + [Sbirths, 0, 0] - Sdie - Sinf;
    I = I + Sinf + Rinf - Idie - Rnew;
    R = R + Rnew - Rinf - Rdie;
    D = D + Sdie + Idie + Rdie;
    N = S + I + R + V;
    
    %% Update tsave, Ssave, Isave, Rsave, Dsave
    tsave(clock) = t; 
    Ssave(clock) = sum(S) ./ sum(N);
    Isave(clock) = sum(I) ./ sum(N);
    Rsave(clock) = sum(R) ./ sum(N);
    Dsave(clock) = sum(D) ./ sum(N);
    Vsave(clock) = sum(V) ./ sum(N);
    % Nsave(clock) = N;

    check = (S+I+R+V) ./ N

    %% Update the plots
    set(hS, 'XData', tsave(1:clock), 'YData', Ssave(1:clock));
    set(hI, 'XData', tsave(1:clock), 'YData', Isave(1:clock));
    set(hR, 'XData', tsave(1:clock), 'YData', Rsave(1:clock));
    set(hD, 'XData', tsave(1:clock), 'YData', Dsave(1:clock));
    set(hV, 'XData', tsave(1:clock), 'YData', Vsave(1:clock));

    % set(hN, 'XData', tsave(1:clock), 'YData', Nsave(1:clock));
    
    drawnow update; % Update the plot
end