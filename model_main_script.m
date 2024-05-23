%=============== Biophysical model for hypoosmotic shock =================%
% Last edited: May 22, 2024
clear, close all

%% Configs

% Generate plots
PLOT_ = 1;

%% Parameters

% Simulation parameters (Euler Method)
t0 = 0;                       % Initial time
tf = 150;                     % Final time
dt = 0.05;                    % Time step
num_steps = (tf - t0) / dt;    % Number of steps

% Time of shock
ts = 50;                   % time when hypoosmotic shock starts

% Flow rate coefficient (permeability)
alpha1 = 0.065;             % Rate coefficient for OM
alpha0 = 0.13;              % Rate coefficient for IM
%alpha0 = 0.008; % aqpZ deletion (low IM permeability)

% Membrane/peptidoglycan stiffness
E_IM_bar = 0.05;           % Stiffness of IM
E_OM_bar = 1.0;            % Stiffness of OM
E_PG_bar = 0.5;            % Stiffness of CW
k0_bar = 5;                % Stiffness of connector molecules (Lpp, ompA, etc.)

% Initial environmental osmolarity
C_out = 0.5;               % Osmotic pressure in the environment

% Hypoosmotic shock osmolarity
C_shock = 0;               % Environmental osmolarity after medium switching

% Intial cellular osmotic pressure
O1i = C_out;               % Initial osmotic pressure in the periplasm
O0i = 1;                   % Initial osmotic pressure in the cytoplasm

% Initial volume of cellular compartments
V0i = 1;                   % Initial volume of cytoplasm
V1i = 0.1;                 % Initial volume of periplasm
Vcwi = V0i;                % Initial "volume" of cell wall
Vtoti = V1i + V0i;         % Initial total volume of the cell

% Rest volume of the cell wall
V_PG_rest = V0i/(1 + (O0i - O1i)/E_PG_bar);    % cell wall rest "volume"

% Lpp softening during OM bulging (including other connector molecules)
g = 0.165;                 % threshold Lpp strain
eps0 = 20;                 % sensitivity of Lpp softening

% MSC parameters
eps1 = 0.1;                % membrane strain threshold
b = 200;                   % sensitivity of MSCs
alpha_MSC = 0.0026;        % osmolyte flow rate coefficient
%alpha_MSC = 0.008;         % mscS/L overexpression

%% Initialization

% Indicator variables
flag_IM = 1; % indicator for physical contact between the IM and CW (whether CW exerts force on IM)

% Timecourse arrays to store results
t = zeros(num_steps, 1);        % time
V1 = zeros(num_steps, 1);       % volume of periplasm
V0 = zeros(num_steps, 1);       % volume of cytoplasm
Vcw = zeros(num_steps, 1);      % volume of cell wall
Vtot = zeros(num_steps, 1);     % total volume of the cell
Vbulge = zeros(num_steps, 1);   % volume of OM bulge
O1 = zeros(num_steps, 1);       % osmotic pressure of periplasm
O0 = zeros(num_steps, 1);       % osmotic pressure of cytoplasm
n0 = zeros(num_steps, 1);       % osmolyte concentration of cytoplasm
n1 = zeros(num_steps, 1);       % osmolyte concentration of periplasm
P1 = zeros(num_steps, 1);       % hydrostatic pressure of the periplasm (periplasm - environment)
P0 = zeros(num_steps, 1);       % hydrostatic pressure of the cytoplasm (cytoplasm - environment)
P_open = zeros(num_steps, 1);    % opening probability of MSCs
e_cyto = zeros(num_steps, 1);      % strain of the inner membrane
Flpp = zeros(num_steps, 1);     % force in connector molecules

% Initial conditions
t(1) = t0;
V1(1) = V1i;
V0(1) = V0i;
Vcw(1) = Vcwi;
Vtot(1) = V1i + V0i;
O1(1) = O1i;
O0(1) = O0i;
P_open(1) = 0;
n0(1) = O0i*V0i; % it is assumed that osmolarity scales with osmolyte concentration
n1(1) = O1i*V1i;

%% Simulation

% Numerical integration using Euler's method
for i = 2:num_steps
    t(i) = t(i-1) + dt;
    
    % Upon hypoosmotic shock (LB to water), decrease environmental osmolarity
    if i==floor(num_steps*(ts/tf))
        C_out = C_shock;
    end   

    % IM strain & MSC opening probability
    e_cyto(i-1) = V0(i-1)/V0i - 1;
    P_open(i-1) = exp(b*(e_cyto(i-1)-eps1))/(1+exp(b*(e_cyto(i-1)-eps1)));
    
    % Forces
    e_peri = (Vtot(i-1) - Vcw(i-1))/V1i - 1; % strain in connector molecules
    Flpp(i) = max([k0_bar * e_peri * exp(-eps0*(e_peri-g))/(1+exp(-eps0*(e_peri-g))),0]); % force in connector molecules; connectors do not resist compression
    Fom = E_OM_bar * (Vtot(i-1)/Vtoti - 1); % OM force
    Fcw = E_PG_bar * (Vcw(i-1)/V_PG_rest - 1); % CW force
    Fim = E_IM_bar * (V0(i-1)/V0i - 1); % IM force

    % Hydrostatic pressure
    P1(i) = Flpp(i) + Fom;
    if flag_IM % pressure is calculated differently depending on whether IM and CW are in contact
        P0(i) = P1(i) + Fim + max([Fcw - Flpp(i), 0]);
    else
        P0(i) = P1(i) + Fim;
    end

    % Water flux rate
    r1 = alpha1 * ((O1(i-1) - C_out) - P1(i));
    r0 = alpha0 * ((O0(i-1) - O1(i-1)) + P1(i) - P0(i));
    
    % Osmolyte flux
    r_MSC = alpha_MSC * P_open(i-1) * (O0(i-1)-O1(i-1)); % rn is positive when: cytoplasm -> periplasm
    n0(i) = n0(i-1) - dt * r_MSC;
    n1(i) = n1(i-1) + dt * r_MSC;
    
    % Water flux & volume changes
    if r_MSC>0 % osmolyte flux out of cytoplasm
        V1(i) = V1(i-1) + dt * r1 - dt * r0 + dt * r_MSC / O0(i-1);
        V0(i) = V0(i-1) + dt * r0 - dt * r_MSC / O0(i-1);
    else % osmolyte flux into cytoplasm
        V1(i) = V1(i-1) + dt * r1 - dt * r0 + dt * r_MSC / O1(i-1);
        V0(i) = V0(i-1) + dt * r0 - dt * r_MSC / O1(i-1);
    end
    Vtot(i) = V1(i) + V0(i);

    % Determine volume of OM bulges
    Vbulge(i) = get_bulge_volume(Vtot(i), V1i, Vcw(i-1), g, eps0);
    
    % Determine cell wall "volume" and whether IM and CW are in contact
    Vb = get_cw_volume(Vtot(i), Vbulge(i), k0_bar, E_PG_bar, V_PG_rest, V1i); % Vb: hypothetical volume of the CW assuming CW and IM are not in contact
    Vcw(i) = max([V0(i), Vb]); % if Vb is larger, take Vb; if Vb is smaller, it means that CW is still in contact with the IM, take V0

    % Osmotic pressure
    O1(i) = n1(i)/V1(i);
    O0(i) = n0(i)/V0(i);
    
    % Update indicator
    if Vb > V0(i)
        flag_IM = 0;
    else
        flag_IM = 1;
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
    plot(t, E_PG_bar * (Vcw./V_PG_rest-1), 'b-', t, E_OM_bar * (Vtot./Vtoti-1), 'r-', t, Flpp, 'g-', t, E_IM_bar * (V0/V0i-1), 'k-');
    xlabel('Time');
    ylabel('Tension');
    legend('Cell wall', 'OM', 'Lpp', 'IM');
    title('Membrane and connector forces');
    xlim([-50,100])
    axis square
    set(gca,'TickDir','out')
    
    % MSC opening
    figure;
    plot(t,e_cyto);
    hold on;
    plot(t,P_open)
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