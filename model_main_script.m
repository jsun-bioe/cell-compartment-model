%=============== Biophysical model for hypoosmotic shock =================%
% Last edited: Jan 11 2024
clear, close all

%% Configs

% Generate plots
PLOT_ = 1;

%% Parameters

% Simulation parameters (Euler Method)
t0 = 0;                       % Initial time
tf = 150;                     % Final time
dt = 0.05;                    % Time step
numSteps = (tf - t0) / dt;    % Number of steps

% Time of shock
ts = 50;                   % time when hypoosmotic shock starts

% Flow rate coefficient (permeability)
alpha1 = 0.065;             % Rate coefficient for OM
alpha0 = 0.13;              % Rate coefficient for IM
%alpha0 = 0.008; % aqpZ deletion (low IM permeability)

% Membrane/peptidoglycan stiffness
kIM = 0.05;                % Stiffness of IM
kOM = 1.0;                   % Stiffness of OM
kCW = 0.5;                 % Stiffness of CW
kLpp = 5;                 % Stiffness of connector molecules (Lpp, ompA, etc.)

% Initial environmental osmolarity
Oout = 0.5;                % Osmotic pressure in the environment

% Hypoosmotic shock osmolarity
Oshock = 0;                % Environmental osmolarity after medium switching

% Intial cellular osmotic pressure
O1i = Oout;                % Initial osmotic pressure in the periplasm
O0i = 1;                   % Initial osmotic pressure in the cytoplasm

% Initial volume of cellular compartments
V0i = 1;                   % Initial volume of cytoplasm
V1i = 0.1;                 % Initial volume of periplasm
Vcwi = V0i;                % Initial "volume" of cell wall
Vtoti = V1i + V0i;         % Initial total volume of the cell

% Rest volume of the cell wall
Vr = V0i/(1 + (O0i - O1i)/kCW);    % cell wall rest "volume"

% Lpp softening during OM bulging (including other connector molecules)
Emax = 0.165;               % threshold Lpp strain
blpp = 20;                 % sensitivity of Lpp softening

% MSC parameters
eps0 = 0.1;                % membrane strain threshold
b = 200;                   % sensitivity of MSCs
alphaM = 0.0026;            % osmolyte flow rate coefficient
%alphaM = 0.008; % mscS/L overexpression

%% Initialization

% Indicator variables
flagIM = 1; % indicator for physical contact between the IM and CW (whether CW exerts force on IM)

% Timecourse arrays to store results
t = zeros(numSteps, 1);        % time
V1 = zeros(numSteps, 1);       % volume of periplasm
V0 = zeros(numSteps, 1);       % volume of cytoplasm
Vcw = zeros(numSteps, 1);      % volume of cell wall
Vtot = zeros(numSteps, 1);     % total volume of the cell
Vbulge = zeros(numSteps, 1);   % volume of OM bulge
O1 = zeros(numSteps, 1);       % osmotic pressure of periplasm
O0 = zeros(numSteps, 1);       % osmotic pressure of cytoplasm
n0 = zeros(numSteps, 1);       % osmolyte concentration of cytoplasm
n1 = zeros(numSteps, 1);       % osmolyte concentration of periplasm
P1 = zeros(numSteps, 1);       % hydrostatic pressure of the periplasm (periplasm - environment)
P0 = zeros(numSteps, 1);       % hydrostatic pressure of the cytoplasm (cytoplasm - environment)
Popen = zeros(numSteps, 1);    % opening probability of MSCs
eps = zeros(numSteps, 1);      % strain of the inner membrane
Flpp = zeros(numSteps, 1);     % force in connector molecules

% Initial conditions
t(1) = t0;
V1(1) = V1i;
V0(1) = V0i;
Vcw(1) = Vcwi;
Vtot(1) = V1i + V0i;
O1(1) = O1i;
O0(1) = O0i;
Popen(1) = 0;
n0(1) = O0i*V0i; % it is assumed that osmolarity scales with osmolyte concentration
n1(1) = O1i*V1i;

%% Simulation

% Numerical integration using Euler's method
for i = 2:numSteps
    t(i) = t(i-1) + dt;
    
    % Upon hypoosmotic shock (LB to water), decrease environmental osmolarity
    if i==floor(numSteps*(ts/tf))
        Oout = Oshock;
    end   

    % IM strain & MSC opening probability
    eps(i-1) = V0(i-1)/V0i - 1;
    Popen(i-1) = exp(b*(eps(i-1)-eps0))/(1+exp(b*(eps(i-1)-eps0)));
    
    % Forces
    elpp = (Vtot(i-1) - Vcw(i-1))/V1i - 1; % strain in connector molecules
    Flpp(i) = max([kLpp * elpp * exp(-blpp*(elpp-Emax))/(1+exp(-blpp*(elpp-Emax))),0]); % force in connector molecules; connectors do not resist compression
    Fom = kOM * (Vtot(i-1)/Vtoti - 1); % OM force
    Fcw = kCW * (Vcw(i-1)/Vr - 1); % CW force
    Fim = kIM * (V0(i-1)/V0i - 1); % IM force

    % Hydrostatic pressure
    P1(i) = Flpp(i) + Fom;
    if flagIM % pressure is calculated differently depending on whether IM and CW are in contact
        P0(i) = P1(i) + Fim + max([Fcw - Flpp(i), 0]);
    else
        P0(i) = P1(i) + Fim;
    end

    % Water flux rate
    r1 = alpha1 * ((O1(i-1) - Oout) - P1(i));
    r0 = alpha0 * ((O0(i-1) - O1(i-1)) + P1(i) - P0(i));
    
    % Osmolyte flux
    rn = alphaM * Popen(i-1) * (O0(i-1)-O1(i-1)); % rn is positive when: cytoplasm -> periplasm
    n0(i) = n0(i-1) - dt * rn;
    n1(i) = n1(i-1) + dt * rn;
    
    % Water flux & volume changes
    if rn>0 % osmolyte flux out of cytoplasm
        V1(i) = V1(i-1) + dt * r1 - dt * r0 + dt * rn / O0(i-1);
        V0(i) = V0(i-1) + dt * r0 - dt * rn / O0(i-1);
    else % osmolyte flux into cytoplasm
        V1(i) = V1(i-1) + dt * r1 - dt * r0 + dt * rn / O1(i-1);
        V0(i) = V0(i-1) + dt * r0 - dt * rn / O1(i-1);
    end
    Vtot(i) = V1(i) + V0(i);

    % Determine volume of OM bulges
    Vbulge(i) = get_bulge_volume(Vtot(i), V1i, Vcw(i-1), Emax, blpp);
    
    % Determine cell wall "volume" and whether IM and CW are in contact
    Vb = get_cw_volume(Vtot(i), Vbulge(i), kLpp, kCW, Vr, V1i); % Vb: hypothetical volume of the CW assuming CW and IM are not in contact
    Vcw(i) = max([V0(i), Vb]); % if Vb is larger, take Vb; if Vb is smaller, it means that CW is still in contact with the IM, take V0

    % Osmotic pressure
    O1(i) = n1(i)/V1(i);
    O0(i) = n0(i)/V0(i);
    
    % Update indicator
    if Vb > V0(i)
        flagIM = 0;
    else
        flagIM = 1;
    end

end

%% Plots
if PLOT_
    t = t-ts;
    
    % Dynamics of compartment volumes
    figure;
    plot(t, V0, 'b-', t, V1, 'r-', t, Vcw, 'g-', t, Vtot, 'k-');
    xlabel('Time');
    ylabel('Volume');
    legend('V0 (Cytoplasm)', 'V1 (Periplasm)', 'Vcw (Cell wall)', 'Vtot (Total volume)');
    title('Dynamics of volumes');
    ylim([0,1.4])
    xlim([-50,100])
    axis square
    set(gca,'TickDir','out')
    
    % Normalized volumes
    figure;
    plot(t, V0./V0(1), 'b-', t, V1./V1(1), 'r-', t, Vtot./Vtot(1), 'k-');
    xlabel('Time');
    ylabel('Nomalized volume');
    legend('Cytoplasm', 'Periplasm', 'Whole cell');
    title('Dynamics of volumes (normalized)');
    ylim([0.75,2.5])
    xlim([-50,100])
    yticks(0.75:0.25:2.5)
    axis square
    set(gca,'TickDir','out')
    
    % Dynamics of osmolytes
    figure;
    plot(t, n0, 'b-', t, n1, 'r-');
    xlabel('Time');
    ylabel('Osmolyte content');
    legend('n0 (Cytoplasm)', 'n1 (Periplasm)');
    title('Dynamics of osmolytes');
    xlim([-50,100])
    axis square
    set(gca,'TickDir','out')
    
    % Dynamics of osmotiv pressure
    figure;
    plot(t, O0, 'b-', t, O1, 'r-', t, O0-O1, 'k-');
    xlabel('Time');
    ylabel('Osmotic pressure');
    legend('Cytoplasm', 'Periplasm', 'Difference');
    title('Dynamics of osmotic pressure');
    xlim([-50,100])
    ylim([0,1])
    axis square
    set(gca,'TickDir','out')
    
    % Forces
    figure;
    plot(t, kCW * (Vcw./Vr-1), 'b-', t, kOM * (Vtot./Vtoti-1), 'r-', t, Flpp, 'g-', t, kIM * (V0/V0i-1), 'k-');
    xlabel('Time');
    ylabel('Tension');
    legend('Cell wall', 'OM', 'Lpp', 'IM');
    title('Membrane and connector forces');
    xlim([-50,100])
    axis square
    set(gca,'TickDir','out')
    
    % MSC opening
    figure;
    plot(t,eps);
    hold on;
    plot(t,Popen)
    hold off;
    title('Dynamics of MSCs: epsilon (blue) and Popen (orange)')
    xlabel('Time');
    ylabel('Eps/Popen');
    xlim([-50,100])
    axis square
    set(gca,'TickDir','out')
    
    % Hydrostatic pressure
    figure, hold on
    plot(t(2:end), P1(2:end), 'r-', t(2:end), P0(2:end), 'b-', t(2:end), P0(2:end)-P1(2:end), 'k-')
    title('Hydrostatic pressure')
    legend('Periplasm', 'Cytoplasm', 'Difference');
    xlabel('Time');
    ylabel('Hydrostatic pressure');
    xlim([-50,100])
    axis square
    set(gca,'TickDir','out')

end

%% Display the final volumes
fprintf('Max volume of cytoplasm (V0): %.2f\n', max(V0));
fprintf('Max volume of periplasm (V1): %.2f\n', max(V1));
fprintf('Final volume of cytoplasm (V0): %.2f\n', V0(end));
fprintf('Final volume of periplasm (V1): %.2f\n', V1(end));