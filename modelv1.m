%% Model Parameters

%% TODO 
%%%% make function 
%%%% realistic parameters deathrate, reinfection, etc
%%%% implement cost (find stats) 

month      = 60*60*24; % Day length (s).
tmax     = month * 40; % Duration of the simulation (s).
clockmax = 100;% Number of time steps.
dt = tmax/clockmax ;% Calculates the duration of each time step.

vDevelopTime = 10 * month;
vDevelopCost = 31 * (10 ^ 9) / vDevelopTime;
vCostPerDose = 20;
qDailyCost = 400;

A           = 0.85/month; 
reprNumb    = 4.3; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7359536/
UB          = A/reprNumb; % https://www.medicalnewstoday.com/articles/how-long-is-a-person-contagious-with-coronavirus
VB          = 1.5*UB;  % https://www.health.com/condition/infectious-diseases/coronavirus/spread-covid-after-vaccine#citation-8:~:text=Some%20research%20has,8

a           = [A, A/2, 0];
b           = [UB, VB, UB];
ra          = 0.5;         % reinfection multiplier

betaH       = (11/1000)/(month * 12);   % birthrate for healthy https://www.cdc.gov/nchs/fastats/births.htm
betaI       = betaH * (1/4); % birthrate for ill 

deltaH      = (798/100000)/(month*12);       % Death rate for healthy individuals https://www.cdc.gov/nchs/products/databriefs/db492.htm#:~:text=Data%20from%20the%20National%20Vital,2021%20to%20798.8%20in%202022.
deltaUI     = 0.083 * UB;       % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9848037/#:~:text=The%20covariate%2Dadjusted%20mortality%20rates%20were%205.1%25%20and%208.3%25%20for%20vaccinated%20and%20unvaccinated%20patients%20hospitalized%20with%20COVID%2D19%2C%20respectively%2C%20in%20the%20whole%20analysis%20sample
deltaVI     = 0.051 * UB;         % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9848037/#:~:text=The%20covariate%2Dadjusted%20mortality%20rates%20were%205.1%25%20and%208.3%25%20for%20vaccinated%20and%20unvaccinated%20patients%20hospitalized%20with%20COVID%2D19%2C%20respectively%2C%20in%20the%20whole%20analysis%20sample
deltaI      = [deltaUI, deltaVI,  deltaUI]; % Death rate for infected individuals

vr          = 0.5/month;       % Vaccination rate % https://usafacts.org/visualizations/covid-vaccine-tracker-states/
qr          = 0.3/month;       % Quarantine rate

initialN = 335 * (10^6); % https://census.gov/quickfacts/fact/table/US/PST045221
intialI = 0.5 * (10^6);

%% Initial Conditions 
% All start with at least one to avoid division by 0 when scaling. 
N =     [initialN,    1,     1] ; % Total population
I =     [intialI,     1,     1] ; % Infected
S =     [N(1)-I(1),   1,     1] ; % Susceptible 
R =     [0,           1,     1] ; % Recovered
D =     [0,           1,     1] ; % Total Deceased

VC = 0;
QC = 0;

%% Initialization for Plotting
tsave = zeros(clockmax, 1);
Ssave = zeros(clockmax, 3);
Isave = zeros(clockmax, 3);
Rsave = zeros(clockmax, 3);
Dsave = zeros(clockmax, 3);
Nsave = zeros(clockmax, 3);
Usave = zeros(clockmax, 3);
Vsave = zeros(clockmax, 3);
Qsave = zeros(clockmax, 3);
NDsave = zeros(clockmax, 1);
NBsave = zeros(clockmax, 1);
VCsave = zeros(clockmax, 1);
QCsave = zeros(clockmax, 1);

%% Create the figure and subplots
figure;
subplot(2,4,1);
hold on
hS1 = plot(tsave(1:clockmax), sum(Ssave(1:clockmax)), 'g', 'LineWidth', 1.5);
hI1 = plot(tsave(1:clockmax), sum(Isave(1:clockmax)), 'r', 'LineWidth', 1.5);
hR1 = plot(tsave(1:clockmax), sum(Rsave(1:clockmax)), 'b', 'LineWidth', 1.5);
hD1 = plot(tsave(1:clockmax), sum(Dsave(1:clockmax)), 'k', 'LineWidth', 1.5);
legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Total Population')

subplot(2,4,2);
hold on
hS2 = plot(tsave(1:clockmax), Ssave(1:clockmax, 1), 'g', 'LineWidth', 1.5);
hI2 = plot(tsave(1:clockmax), Isave(1:clockmax, 1), 'r', 'LineWidth', 1.5);
hR2 = plot(tsave(1:clockmax), Rsave(1:clockmax, 1), 'b', 'LineWidth', 1.5);
hD2 = plot(tsave(1:clockmax), Dsave(1:clockmax, 1), 'k', 'LineWidth', 1.5);
legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Not vaccinated')

subplot(2,4,3);
hold on
hS3 = plot(tsave(1:clockmax), Ssave(1:clockmax, 2), 'g', 'LineWidth', 1.5);
hI3 = plot(tsave(1:clockmax), Isave(1:clockmax, 2), 'r', 'LineWidth', 1.5);
hR3 = plot(tsave(1:clockmax), Rsave(1:clockmax, 2), 'b', 'LineWidth', 1.5);
hD3 = plot(tsave(1:clockmax), Dsave(1:clockmax, 2), 'k', 'LineWidth', 1.5);
legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Vaccinated')

subplot(2,4,4);
hold on
hS4 = plot(tsave(1:clockmax), Ssave(1:clockmax, 3), 'g', 'LineWidth', 1.5);
hI4 = plot(tsave(1:clockmax), Isave(1:clockmax, 3), 'r', 'LineWidth', 1.5);
hR4 = plot(tsave(1:clockmax), Rsave(1:clockmax, 3), 'b', 'LineWidth', 1.5);
hD4 = plot(tsave(1:clockmax), Dsave(1:clockmax, 3), 'k', 'LineWidth', 1.5);
legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Quarantined')

subplot(2,4,5)
hold on
hU  = plot(tsave(1:clockmax), Usave(1:clockmax), 'g', 'LineWidth', 1.5);
hV  = plot(tsave(1:clockmax), Vsave(1:clockmax), 'r', 'LineWidth', 1.5);
hQ  = plot(tsave(1:clockmax), Qsave(1:clockmax), 'k', 'LineWidth', 1.5);
legend({'U','V','Q'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('States')

subplot(2,4,6);
hold on
hN  = plot(tsave(1:clockmax), sum(Nsave(1:clockmax), 2), 'c', 'LineWidth', 1.5);
hND = plot(tsave(1:clockmax), NDsave(1:clockmax), 'k', 'LineWidth', 1.5);
hNB = plot(tsave(1:clockmax), NBsave(1:clockmax), 'g', 'LineWidth', 1.5);
expectedSize = sum(N) * (1 + betaH - deltaH)^tmax;
legend({'N','D','B'},'Location','northeast')
axis([0, tmax/month, 0, expectedSize]);
title('Total Population')

subplot(2,4,7);
hold on
hTC  = plot(tsave(1:clockmax), QCsave(1:clockmax), 'k', 'LineWidth', 1.5);
hVC = plot(tsave(1:clockmax), VCsave(1:clockmax), 'r', 'LineWidth', 1.5);
hQC = plot(tsave(1:clockmax), QCsave(1:clockmax), 'b', 'LineWidth', 1.5);
expectedSize = sum(N) * (1 + betaH - deltaH)^tmax;
legend({'TC','VC','QC'},'Location','northeast')
axis([0, tmax/month, 0, 2.5*10^11]);
title('Costs')

drawnow;

births = 0;

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

    %% Quarantines
    S1toS3 = S(1) * qr * dt;
    R1toR3 = R(1) * qr * dt;
    I1toI3 = I(1) * qr*2 * dt;
    
    S(1) = S(1) - S1toS3;
    S(3) = S(3) + S1toS3;
    
    I(1) = I(1) - I1toI3;
    I(3) = I(3) + I1toI3;
    
    R(1) = R(1) - R1toR3;
    R(3) = R(3) + R1toR3;

    QC = QC + qDailyCost * (newS3 + newI3 + newR3);

    %% Vaccinations
    if t > vDevelopTime
        R1toR2 = R(1) * vr * dt;
        R3toR2 = R(3) * vr * dt;
        S1toS2 = S(1) * vr * dt;
        S3toS2 = S(3) * vr * dt;
        VC = VC + vCostPerDose * (newR2 + newS2);
    else
        R1toR2 = 0;
        R3toR2 = 0;
        S1toS2 = 0;
        S3toS2 = 0;
        
        VC = t * vDevelopCost;
    end

    S(1) = S(1) - R1toR2;
    S(2) = S(2) + R1toR2 + R3toR2;
    S(3) = S(3) - R3toR2;
    
    R(1) = R(1) - R1toR2;
    R(2) = R(2) + R1toR2 + R3toR2;
    R(3) = R(3) - R3toR2;
    
    %% Results
    N = S + R + I;
    births = births + Sbirths;

    if t < vDevelopTime
        N(2) = 1;
        S(2) = 0;
        I(2) = 0;
        R(2) = 0;
        D(2) = 0;
    end
        

    % Update tsave, Ssave, Isave, Rsave, Dsave
    tsave(clock) = t / month; 
    Nsave(clock, :) = N;
    Ssave(clock, :) = S;
    Isave(clock, :) = I;
    Rsave(clock, :) = R;
    Dsave(clock, :) = D;
    Usave(clock, :) = N(1) / sum(N);
    Vsave(clock, :) = N(2) / sum(N);
    Qsave(clock, :) = N(3) / sum(N);
    NDsave(clock) = sum(D);
    NBsave(clock) = births;

    VCsave(clock) = VC;
    QCsave(clock) = QC;
    
    % Update the plots in the first subplot
    subplot(2,4,1);
    scale = sum(Nsave(1:clock, :), 2);
    set(hS1, 'XData', tsave(1:clock), 'YData', sum(Ssave(1:clock, :) ./ scale, 2));
    set(hI1, 'XData', tsave(1:clock), 'YData', sum(Isave(1:clock, :) ./ scale, 2));
    set(hR1, 'XData', tsave(1:clock), 'YData', sum(Rsave(1:clock, :) ./ scale, 2));
    set(hD1, 'XData', tsave(1:clock), 'YData', sum(Dsave(1:clock, :) ./ scale, 2));


    % Update the plots in the second subplot
    subplot(2,4,2);
    scale = Nsave(1:clock, 1);
    set(hS2, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 1) ./ scale);
    set(hI2, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 1) ./ scale);
    set(hR2, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 1) ./ scale);
    set(hD2, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 1) ./ scale);

     % Update the plots in the first subplot
    subplot(2,4,3);
    scale = Nsave(1:clock, 2);
    set(hS3, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 2) ./ scale);
    set(hI3, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 2) ./ scale);
    set(hR3, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 2) ./ scale);
    set(hD3, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 2) ./ scale);
    
    % Update the plots in the second subplot
    subplot(2,4,4);
    scale = Nsave(1:clock, 3);
    set(hS4, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 3) ./ scale);
    set(hI4, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 3) ./ scale);
    set(hR4, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 3) ./ scale);
    set(hD4, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 3) ./ scale);
    
    subplot(2,4,5);
    set(hU, 'XData',  tsave(1:clock), 'YData', Usave(1:clock));
    set(hV, 'XData',  tsave(1:clock), 'YData', Vsave(1:clock));
    set(hQ, 'XData',  tsave(1:clock), 'YData', Qsave(1:clock));

    subplot(2,4,6);
    set(hN, 'XData',  tsave(1:clock), 'YData', sum(Nsave(1:clock, :), 2));
    set(hND, 'XData',  tsave(1:clock), 'YData', NDsave(1:clock));
    set(hNB, 'XData',  tsave(1:clock), 'YData', NBsave(1:clock));

    subplot(2,4,7);
    set(hTC, 'XData',  tsave(1:clock), 'YData', QCsave(1:clock) + VCsave(1:clock));
    set(hQC, 'XData',  tsave(1:clock), 'YData', QCsave(1:clock));
    set(hVC, 'XData',  tsave(1:clock), 'YData', VCsave(1:clock));

    drawnow;
    
end
