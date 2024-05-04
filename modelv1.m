%% Model Parameters

%% TODO 
%%%% make function 
%%%% realistic parameters deathrate, reinfection, etc
%%%% implement cost (find stats) 

imax = 0; 

m = 5;
day      = 24; % day length (h).
month    = day*30; % Month length (h).
tmax     = month * 40; % Duration of the simulation (h).
clockmax = 100;% Number of time steps.
dt = tmax/clockmax ;% Calculates the duration of each time step.

vDevelopTime = 10 * month;
vDevelopCost = 31 * (10 ^ 9) / vDevelopTime;
vCostPerDose = 20;
qDailyCostI = 70/day;
qDailyCostH = 50/day;

A           = 0.85/month; 
reprNumb    = 4.3; % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7359536/
UB          = A/reprNumb; % https://www.medicalnewstoday.com/articles/how-long-is-a-person-contagious-with-coronavirus
VB          = 1.5*UB;  % https://www.health.com/condition/infectious-diseases/coronavirus/spread-covid-after-vaccine#citation-8:~:text=Some%20research%20has,8

a           = [A, A/2, 0];
b           = [UB, VB, UB];
ra          = 0.5;         % reinfection multiplier

betaH       = (11/1000)/(month * 12);   % birthrate for healthy https://www.cdc.gov/nchs/fastats/births.htm
betaI       = betaH * (1/4);            % birthrate for ill 

deltaH      = (798/100000)/(month*12);       % Death rate for healthy individuals https://www.cdc.gov/nchs/products/databriefs/db492.htm#:~:text=Data%20from%20the%20National%20Vital,2021%20to%20798.8%20in%202022.
deltaUI     = 0.083 / month;         % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9848037/#:~:text=The%20covariate%2Dadjusted%20mortality%20rates%20were%205.1%25%20and%208.3%25%20for%20vaccinated%20and%20unvaccinated%20patients%20hospitalized%20with%20COVID%2D19%2C%20respectively%2C%20in%20the%20whole%20analysis%20sample
deltaVI     = 0.051 / month;         % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9848037/#:~:text=The%20covariate%2Dadjusted%20mortality%20rates%20were%205.1%25%20and%208.3%25%20for%20vaccinated%20and%20unvaccinated%20patients%20hospitalized%20with%20COVID%2D19%2C%20respectively%2C%20in%20the%20whole%20analysis%20sample
deltaI      = [deltaUI, deltaVI, deltaUI]; % Death rate for infected individuals

vr          = 0.1 /month;       % Vaccination rate % https://usafacts.org/visualizations/covid-vaccine-tracker-states/
qr          = .01/month;        % Quarantine rate

intialI = 0.5 * (10^6);
initialN = 335 * (10^6) + intialI; % https://census.gov/quickfacts/fact/table/US/PST045221

%% Initial Conditions 
% All start with at least one to avoid division by 0 when scaling. 
N =     [initialN,    1,     1] ; % Total population
I =     [intialI,     1,     1] ; % Infected
S =     [N(1)-I(1),   1,     1] ; % Susceptible 
R =     [1,           1,     1] ; % Recovered
D =     [1,           1,     1] ; % Total Deceased

VC = 0;
QC = 0;

CVC = 0;
CQC = 0;

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
CVCsave = zeros(clockmax, 1);
CQCsave = zeros(clockmax, 1);



%% Create the figure and subplots
figure('NumberTitle','off', 'Name', [' vr ' num2str(vr*month) 'qr' num2str(qr*month)]);

subplot(2,m,1);
hold on
hS1 = plot(tsave(1:clockmax), sum(Ssave(1:clockmax)), 'g', 'LineWidth', 1.5);
hI1 = plot(tsave(1:clockmax), sum(Isave(1:clockmax)), 'r', 'LineWidth', 1.5);
hR1 = plot(tsave(1:clockmax), sum(Rsave(1:clockmax)), 'b', 'LineWidth', 1.5);
hD1 = plot(tsave(1:clockmax), sum(Dsave(1:clockmax)), 'k', 'LineWidth', 1.5);
legend({'S','I','R', 'D'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Total Population SIRD')

subplot(2,m,2);
hold on
hS2 = plot(tsave(1:clockmax), Ssave(1:clockmax, 1), 'g', 'LineWidth', 1.5);
hI2 = plot(tsave(1:clockmax), Isave(1:clockmax, 1), 'r', 'LineWidth', 1.5);
hR2 = plot(tsave(1:clockmax), Rsave(1:clockmax, 1), 'b', 'LineWidth', 1.5);
%hD2 = plot(tsave(1:clockmax), Dsave(1:clockmax, 1), 'k', 'LineWidth', 1.5);
legend({'S','I','R'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Not vaccinated')

subplot(2,m,3);
hold on
hS3 = plot(tsave(1:clockmax), Ssave(1:clockmax, 2), 'g', 'LineWidth', 1.5);
hI3 = plot(tsave(1:clockmax), Isave(1:clockmax, 2), 'r', 'LineWidth', 1.5);
hR3 = plot(tsave(1:clockmax), Rsave(1:clockmax, 2), 'b', 'LineWidth', 1.5);
%hD3 = plot(tsave(1:clockmax), Dsave(1:clockmax, 2), 'k', 'LineWidth', 1.5);
legend({'S','I','R'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Vaccinated')

subplot(2,m,4);
hold on
hS4 = plot(tsave(1:clockmax), Ssave(1:clockmax, 3), 'g', 'LineWidth', 1.5);
hI4 = plot(tsave(1:clockmax), Isave(1:clockmax, 3), 'r', 'LineWidth', 1.5);
hR4 = plot(tsave(1:clockmax), Rsave(1:clockmax, 3), 'b', 'LineWidth', 1.5);
%hD4 = plot(tsave(1:clockmax), Dsave(1:clockmax, 3), 'k', 'LineWidth', 1.5);
legend({'S','I','R'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('Quarantined')

subplot(2,m,5)
hold on
hU  = plot(tsave(1:clockmax), Usave(1:clockmax), 'g', 'LineWidth', 1.5);
hV  = plot(tsave(1:clockmax), Vsave(1:clockmax), 'r', 'LineWidth', 1.5);
hQ  = plot(tsave(1:clockmax), Qsave(1:clockmax), 'k', 'LineWidth', 1.5);
legend({'U','V','Q'},'Location','northeast')
axis([0, tmax/month, 0, 1.02])
title('States')

subplot(2,m,6);
hold on
hN  = plot(tsave(1:clockmax), sum(Nsave(1:clockmax), 2), 'c', 'LineWidth', 1.5);
hND = plot(tsave(1:clockmax), NDsave(1:clockmax), 'k', 'LineWidth', 1.5);
hNB = plot(tsave(1:clockmax), NBsave(1:clockmax), 'g', 'LineWidth', 1.5);
expectedSize = sum(N) * (1 + betaH - deltaH)^tmax;
legend({'N','D','B'},'Location','northeast')
axis([0, tmax/month, 0, expectedSize]);
title('Total Population Size Change')

subplot(2,m,8);
hold on
hTC  = plot(tsave(1:clockmax), QCsave(1:clockmax), 'k', 'LineWidth', 1.5);
hVC = plot(tsave(1:clockmax), VCsave(1:clockmax), 'r', 'LineWidth', 1.5);
hQC = plot(tsave(1:clockmax), QCsave(1:clockmax), 'b', 'LineWidth', 1.5);
expectedSize = sum(N) * (1 + betaH - deltaH)^tmax;
legend({'TC','VC','QC'},'Location','northeast')
axis([0, tmax/month, 0, 5*10^10]);
title('Cost at Time')

subplot(2,m,7);
hold on
hDT = plot(tsave(1:clockmax), sum(Dsave(1:clockmax), 2), 'k', 'LineWidth', 1.5);
hDU = plot(tsave(1:clockmax), Dsave(1:clockmax, 1), 'r', 'LineWidth', 1.5);
hDV = plot(tsave(1:clockmax), Dsave(1:clockmax, 2), 'g', 'LineWidth', 1.5);
hDQ = plot(tsave(1:clockmax), Dsave(1:clockmax, 3), 'b', 'LineWidth', 1.5);
legend({'TD', 'UD','VD','QD'},'Location','northeast')
axis([0, tmax/month, 0, sum(D)]);
title('Deaths')

subplot(2,m,9);
 hold on
 hCTC = plot(tsave(1:clockmax), CQCsave(1:clockmax), 'k', 'LineWidth', 1.5);
 hCVC = plot(tsave(1:clockmax), CVCsave(1:clockmax), 'r', 'LineWidth', 1.5);
 hCQC = plot(tsave(1:clockmax), CQCsave(1:clockmax), 'b', 'LineWidth', 1.5);
 legend({'CTC','CVC','CQC'},'Location','northeast')
 axis([0, tmax/month, 0, 9*10^11]);
 title('Cumulative Costs')




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
    if (I>imax)
        imax = I;
    end
    infprop = sum(I)/sum(N);
    %% Quarantines
    S1toS3 = S(1) * qr * infprop * dt;
    R1toR3 = R(1) * qr * infprop * dt;
    I1toI3 = I(1) * qr*2 * infprop * dt;
    
    S(1) = S(1) - S1toS3;
    S(3) = S(3) + S1toS3;
    
    I(1) = I(1) - I1toI3;
    I(3) = I(3) + I1toI3;
    
    R(1) = R(1) - R1toR3;
    R(3) = R(3) + R1toR3;

    %% Vaccinations
    if t > vDevelopTime
        R1toR2 = R(1) * vr * dt;
        R3toR2 = R(3) * vr * dt;
        S1toS2 = S(1) * vr * dt;
        S3toS2 = S(3) * vr * dt;
        VC = vCostPerDose * (R1toR2 + R3toR2 + S1toS2 + S3toS2);
        CVC = CVC + VC;
    else
        R1toR2 = 0;
        R3toR2 = 0;
        S1toS2 = 0;
        S3toS2 = 0;
        
        VC = t * vDevelopCost;
        CVC = VC;
        % totalDevCost = VC;
    end

    S(1) = S(1) - S1toS2;
    S(2) = S(2) + S1toS2 + S3toS2;
    S(3) = S(3) - S3toS2;
    
    R(1) = R(1) - R1toR2;
    R(2) = R(2) + R1toR2 + R3toR2;
    R(3) = R(3) - R3toR2;
    
    %% Results
    N = S + R + I;
    births = births + Sbirths;

    QC =  I(3) * qDailyCostI * dt + (S(3) + R(3)) * qDailyCostH *dt;
    CQC = CQC + QC;
    if t < vDevelopTime
        N(2) = 1;
        S(2) = 0;
        I(2) = 0;
        R(2) = 0;
        D(2) = 0;
    end

    if qr == 0
        N(3) = 1;
        S(3) = 0;
        I(3) = 0;
        R(3) = 0;
        D(3) = 0;
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
    CVCsave(clock) = CVC;
    CQCsave(clock) = CQC;
    
    % Update the plots in the first subplot
    subplot(2,m,1);
    scale = sum(Nsave(1:clock, :), 2);
    set(hS1, 'XData', tsave(1:clock), 'YData', sum(Ssave(1:clock, :) ./ scale, 2));
    set(hI1, 'XData', tsave(1:clock), 'YData', sum(Isave(1:clock, :) ./ scale, 2));
    set(hR1, 'XData', tsave(1:clock), 'YData', sum(Rsave(1:clock, :) ./ scale, 2));
    set(hD1, 'XData', tsave(1:clock), 'YData', sum(Dsave(1:clock, :) ./ scale, 2));


    % Update the plots in the second subplot
    subplot(2,m,2);
    scale = Nsave(1:clock, 1);
    set(hS2, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 1) ./ scale);
    set(hI2, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 1) ./ scale);
    set(hR2, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 1) ./ scale);
    %set(hD2, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 1) ./ scale);

     % Update the plots in the first subplot
    subplot(2,m,3);
    scale = Nsave(1:clock, 2);
    set(hS3, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 2) ./ scale);
    set(hI3, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 2) ./ scale);
    set(hR3, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 2) ./ scale);
    %set(hD3, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 2) ./ scale);
    
    % Update the plots in the second subplot
    subplot(2,m,4);
    scale = Nsave(1:clock, 3);
    set(hS4, 'XData', tsave(1:clock), 'YData', Ssave(1:clock, 3) ./ scale);
    set(hI4, 'XData', tsave(1:clock), 'YData', Isave(1:clock, 3) ./ scale);
    set(hR4, 'XData', tsave(1:clock), 'YData', Rsave(1:clock, 3) ./ scale);
    %set(hD4, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 3) ./ scale);
    
    subplot(2,m,5);
    set(hU, 'XData',  tsave(1:clock), 'YData', Usave(1:clock));
    set(hV, 'XData',  tsave(1:clock), 'YData', Vsave(1:clock));
    set(hQ, 'XData',  tsave(1:clock), 'YData', Qsave(1:clock));

    subplot(2,m,6);
    set(hN, 'XData',  tsave(1:clock), 'YData', sum(Nsave(1:clock, :), 2));
    set(hND, 'XData',  tsave(1:clock), 'YData', NDsave(1:clock));
    set(hNB, 'XData',  tsave(1:clock), 'YData', NBsave(1:clock));

    subplot(2,m,8);
    set(hTC, 'XData',  tsave(1:clock), 'YData', QCsave(1:clock) + VCsave(1:clock));
    set(hQC, 'XData',  tsave(1:clock), 'YData', QCsave(1:clock));
    set(hVC, 'XData',  tsave(1:clock), 'YData', VCsave(1:clock));

    subplot(2,m,7);
    axis([0, tmax/month, 0, sum(D)*1.05]);
    set(hDT, 'XData',  tsave(1:clock), 'YData', NDsave(1:clock));
    set(hDU, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 1));
    set(hDV, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 3));
    set(hDQ, 'XData', tsave(1:clock), 'YData', Dsave(1:clock, 2));

    subplot(2,m,9);
    set(hCTC, 'XData',  tsave(1:clock), 'YData', CQCsave(1:clock) + CVCsave(1:clock));
    set(hCQC, 'XData',  tsave(1:clock), 'YData', CQCsave(1:clock));
    set(hCVC, 'XData',  tsave(1:clock), 'YData', CVCsave(1:clock));

    drawnow;
    
end