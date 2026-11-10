clc; clear all; close all;

%% ==================== CONFIGURABLE PARAMETERS ====================
N_segments = 2;

% Segment 1 (Horizontal) and Segment 2 (Vertical)
L     = [0.5, 0.5];       % Length of each segment (m)
E     = [200e9, 200e9];   % Young's Modulus (Pa)
Width = [20e-3, 20e-3];   % Cross-section width (m)
Thk   = [2e-3, 2e-3];     % Cross-section thickness (m)
rho   = [7850, 7850];     % Density (kg/m^3)

% Frequency sweep range for finding natural frequencies (rad/s)
w_test = 1:1:1000;
%% =================================================================

% Derived Geometric Quantities
I = (1/12) .* Width .* Thk.^3; % Second moment of area (m^4)
A = Width .* Thk;              % Cross-sectional area (m^2)

%% TASK 2: Point Transfer Matrix for 90-degree Corner
R = [ 0,  1,  0,  0,  0,  0;  % u2 = w1
     -1,  0,  0,  0,  0,  0;  % w2 = -u1
      0,  0,  1,  0,  0,  0;  % theta2 = theta1
      0,  0,  0,  0,  1,  0;  % N2 = V1
      0,  0,  0, -1,  0,  0;  % V2 = -N1
      0,  0,  0,  0,  0,  1]; % M2 = M1

%% Find the natural frequencies
D = zeros(1, length(w_test));

for i = 1:length(w_test)
    w = w_test(i);
    GTM = eye(6);

    for s = 1:N_segments
        B = zeros(6,6);
        B(1,4) = 1/(E(s)*A(s));
        B(2,3) = 1;
        B(3,6) = 1/(E(s)*I(s));
        B(4,1) = -rho(s)*A(s)*w^2;
        B(5,2) = -rho(s)*A(s)*w^2;
        B(6,5) = -1;

        FTM = expm(B * L(s));

        if s == 2
            GTM = R * GTM;
        end

        GTM = FTM * GTM;
    end

    U_sub = GTM(4:6, 4:6);
    D(i) = det(U_sub);
end

% Get Natural Frequencies (Zero Crossings)
sign_changes = find(diff(sign(D)) ~= 0);
num_modes = min(3, length(sign_changes));
omega_n = w_test(sign_changes(1:num_modes));

disp('Natural Frequencies for L-Shaped Fixed-Free Beam:');
fprintf('  (rad/s): %s\n', mat2str(omega_n, 5));
fprintf('  (Hz)   : %s\n', mat2str(omega_n * 0.5 / pi, 5));

%% Solve for mode shapes from GTM and plot in SEPARATE FIGURES
colors = ['r', 'g', 'b'];

for i = 1:num_modes
    w = omega_n(i);

    % Re-evaluate individual FTMs at the exact natural frequency
    B1 = zeros(6,6); B1(1,4)=1/(E(1)*A(1)); B1(2,3)=1; B1(3,6)=1/(E(1)*I(1)); B1(4,1)=-rho(1)*A(1)*w^2; B1(5,2)=-rho(1)*A(1)*w^2; B1(6,5)=-1;
    FTM1 = expm(B1 * L(1));

    B2 = zeros(6,6); B2(1,4)=1/(E(2)*A(2)); B2(2,3)=1; B2(3,6)=1/(E(2)*I(2)); B2(4,1)=-rho(2)*A(2)*w^2; B2(5,2)=-rho(2)*A(2)*w^2; B2(6,5)=-1;
    FTM2 = expm(B2 * L(2));

    GTM = FTM2 * R * FTM1;

    % Find starting state vector (Z0) at the wall
    U_sub = GTM(4:6, 4:6);
    forces = -U_sub(1:2, 1:2) \ U_sub(1:2, 3) * 1;
    N0 = forces(1); V0 = forces(2); M0 = 1;
    Z = [0; 0; 0; N0; V0; M0];

    % Walk down Segment 1 (Horizontal)
    x1 = linspace(0, L(1), 30);
    dX1 = zeros(size(x1)); dY1 = zeros(size(x1));
    for k = 1:length(x1)
        U_x = expm(B1 * x1(k));
        Z_x = U_x * Z;
        dX1(k) = Z_x(1);
        dY1(k) = Z_x(2);
    end

    % Pass state through the corner
    Z = FTM1 * Z;
    Z = R * Z;

    % Walk down Segment 2 (Vertical)
    x2 = linspace(0, L(2), 30);
    dX2 = zeros(size(x2)); dY2 = zeros(size(x2));
    for k = 1:length(x2)
        U_x = expm(B2 * x2(k));
        Z_x = U_x * Z;
        dX2(k) = -Z_x(2);
        dY2(k) = Z_x(1);
    end

    % Scale deflections for visualization
    scale = 0.25 * L(1) / max(max(abs([dX1, dX2, dY1, dY2])));

    % --- CREATE A NEW SEPARATE FIGURE FOR EACH MODE ---
    figure('Color', 'w', 'Name', ['Mode Shape ', num2str(i)]);
    hold on;

    % Plot undeformed reference line
    plot([0, L(1), L(1)], [0, 0, L(2)], 'k--', 'linewidth', 1.5, 'DisplayName', 'Undeformed');

    % Plot deformed mode shape
    plot(x1 + scale*dX1, zeros(size(x1)) + scale*dY1, colors(i), 'linewidth', 2, 'DisplayName', ['Mode ', num2str(i)]);
    plot(L(1)*ones(size(x2)) + scale*dX2, x2 + scale*dY2, colors(i), 'linewidth', 2, 'HandleVisibility', 'off');

    % Formatting options for clear comparison
    title(sprintf('TMM Mode Shape %d  (f = %.2f Hz, \\omega = %d rad/s)', i, omega_n(i)*0.5/pi, omega_n(i)));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    axis equal; grid on; box on; legend('Location', 'best');
end
