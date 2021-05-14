
function NACA_solver()

    clc;

    f = 2;
    p = 4;
    n = 100;
    xh = -1;
    eta = 10;
    u_inf = 1;
    alpha = 0;
    x_ref = 0.25;
    c = 1;
    
    %discussion();
    
    %checkFlapEfficiency(f, p, [0.9]);

    %[x, z] = getAirfoilChamberPoints(f, p, n, xh, eta);
    
    %gamma = solveAirfoil(x, z, alpha, u_inf);
    
    %cl = getLiftCoeficient(gamma, u_inf, c)
    %cm = getMomentCoeficient(gamma, u_inf, alpha, x, x_ref)
    
    [cl_slope, alpha_0, cm0] = getAirfoilCharacteristics(f, p, xh, eta)
    
    %testConvergence(f, p, alpha, xh, eta, u_inf, x_ref, 1.3134, -0.0515);
   
    %cp = getPressureCoeficients(gamma, u_inf, x, z);
    
    %hold on;
    %plot(x, [cp, 0], x, z);
    %legend('pressure', 'profile');
    %grid on;
    %hold off;

end


% Given a 4 digit NACA airfoil described by the parameters
%
%   f: maximum chamber in percent of chord (min 0 max 9)
%   p: maximum chamber position in tenths of chord (min 0 max 9)
%   n: number of pannels used (min 2)
%   xh: flap hinge location in tenths of chord (min 0 max 10) use the
%       a value of 10 or negative ones if no flap is desired
%   eta: flap deflection angle in degree (min 0 max 90)
%
% It returns the positions of the union points between pannels (including
% position of the first and last points)
function [x, z] = getAirfoilChamberPoints(f, p, n, xh, eta)

    f = f/100;
    p = p/10;
    xh = xh/10;
    eta = eta*pi/180;
    
    
    theta = linspace(0, pi, n+1);
    x = 1/2*(1 - cos(theta));
    le_part = (x < p);   % 0 <= x < p
    te_part = (x >= p);  % p <= x <= 1
    z = zeros(1, n+1);
    z(le_part) = f/p^2 * (2*p*x(le_part) - x(le_part).^2);
    z(te_part) = f/(1-p)^2 * (1 - 2*p + 2*p*x(te_part) - x(te_part).^2);
    
    % if flap is desired
    if xh < 1 && xh > 0
        flap_part = (x >= xh);
        z(flap_part) = z(flap_part) - tan(eta)*(x(flap_part) - xh);
    end

end


% Given an airfoil described by the parameters
%
%   x: vector containing the x position of the panel unions that compose
%      the chord of the profile (all values must be between 0 and 1)
%   z: vector containing the z position of the panel unions that compose
%      the chord of the profile (x and z must be the same size)
%   alpha: angle of attack in degrees (recomended min -10 max 10)
%   u_inf: free airstream speed
%
% it returns the values of the vortex intensity for each panel
function gamma = solveAirfoil(x, z, alpha, u_inf)

    % get number of panels used
    n = size(x, 2) - 1;
    alpha = alpha*pi/180;
    
    % compute normal vectors to the pannels
    delta_x = x(2:n+1) - x(1:n);
    delta_z = z(2:n+1) - z(1:n);
    c = sqrt(delta_x.^2 + delta_z.^2);
    nx = -delta_z./c;
    nz = delta_x./c;
    
    % get the positions of the points where the velocity is
    % evaluated (x_v 3c/4) and where the vortex is located (x_ 0 c/4)
    x_v = (delta_x)*3/4 + x(1:n);
    z_v = (delta_z)*3/4 + z(1:n);
    x_0 = (delta_x)*1/4 + x(1:n);
    z_0 = (delta_z)*1/4 + z(1:n);
    
    % prepare system of equations
    lhs = zeros(n);     % right hand side of the equation
    rhs = zeros(n, 1);  % left hand side of the equation
    for i = 1:n
        for j = 1:n
            r_sqrd = (x_v(i) - x_0(j))^2 + (z_v(i) - z_0(j))^2;
            u = 1/(2*pi) * (z_v(i) - z_0(j))/r_sqrd;
            w = -1/(2*pi) * (x_v(i) - x_0(j))/r_sqrd;
            
            lhs(i,j) = dot([u w], [nx(i) nz(i)]);
        end
        
        rhs(i) = -u_inf * dot([cos(alpha) sin(alpha)], [nx(i) nz(i)]);
    end
        
    % solve system of equations and get the values of gamma
    gamma = linsolve(lhs, rhs);

end


% Given the following values of a wing profile
%
%   gamma: vortex density over the chord of the profile
%   u_inf: free airstream speed
%   c: size of the chord
%
% it returns the value of the lift coeficient
function cl = getLiftCoeficient(gamma, u_inf, c)

    cl = 2*sum(gamma)/(u_inf*c);

end


% Given the following values of a wing profile
%
%   gamma: vortex density over the chord of the profile
%   u_inf: free airstream speed
%   alpha: angle of attack of the profile in degrees
%   x: x coordinate position for each panel union (including first and
%      last points)
%   x_ref: reference point for the moment of inertia
%
% returns the value of the moment of inertia
function cm = getMomentCoeficient(gamma, u_inf, alpha, x, x_ref)

    alpha = alpha*pi/180;
    n = size(x, 2) - 1;   % number of panels used
    x_vortex = (x(2:n+1) - x(1:n))*1/4 + x(1:n);
    
    cm = -2*sum(gamma'.*(x_vortex - x_ref)*cos(alpha))/(u_inf*x(n+1)^2);

end


% Given the following values of a wing profile
%
%   gamma: vortex density over the chord of the profile
%   u_inf: free airstream speed
%   alpha: angle of attack of the profile in degrees
%   x: x coordinate position for each panel union (including first and
%      last points)
%   z: z coordinate position for each panel union (including first and
%      last points)
%
% returns the pressure coeficient (Cp) for each panel
function cp = getPressureCoeficients(gamma, u_inf, x, z)

    n = size(x, 2) - 1;   % number of panels used
    delta_x = x(2:n+1) - x(1:n);
    delta_z = z(2:n+1) - z(1:n);
    c = sqrt(delta_x.^2 + delta_z.^2);
    
    cp = 2*gamma'./(u_inf*c);

end


function testConvergence(f, p, alpha, xh, eta, u_inf, x_ref, real_cl, real_cm)

    n = [2 4 8 16 32 64 128 256 512 1024];

    cl = zeros(size(n));
    cm = zeros(size(n));
    
    for i = 1:size(n, 2)
        [x, z] = getAirfoilChamberPoints(f, p, n(i), xh, eta);
        gamma = solveAirfoil(x, z, alpha, u_inf);
        cl(i) = getLiftCoeficient(gamma, u_inf, 1);
        cm(i) = getMomentCoeficient(gamma, u_inf, alpha, x, x_ref);
    end
    
    error_cl = cl - real_cl;
    error_cm = cm - real_cm;
    
    log_n = log2(n);
    
    
    subplot(1, 2, 1);
    yyaxis left;
    plot(log_n, cl);
    ylabel('Cl value');
    yyaxis right;
    plot(log_n, error_cl);
    title('Cl value and error');
    xlabel('log2 number of panels');
    ylabel('error');
    lgd = legend('Cl values', 'error');
    lgd.Location = 'southeast';
    grid on;
    
    subplot(1, 2, 2);
    yyaxis left;
    plot(log_n, cm);
    ylabel('Cm value');
    yyaxis right;
    plot(log_n, error_cm);
    title('Cm value and error');
    xlabel('log2 number of panels');
    ylabel('error');
    legend('Cm values', 'error');
    grid on;

end


function [cl_slope, alpha_0, cm0] = getAirfoilCharacteristics(f, p, xh, eta)

    n = 256; % this value should be big enough
    u_inf = 1;
    
    cl = zeros(1, 10);
    cm = zeros(1, 10);
    
    [x, z] = getAirfoilChamberPoints(f, p, n, xh, eta);
    
    % get a bounch of points for the airfoil Cl
    alpha = linspace(-10, 10, 10);
    
    for i = 1:10
        gamma = solveAirfoil(x, z, alpha(i), u_inf);
        cl(i) = getLiftCoeficient(gamma, u_inf, 1);
        cm(i) = getMomentCoeficient(gamma, u_inf, alpha(i), x, 0.25);
    end
    
    cm0 = mean(cm);
    pol_fit = polyfit(alpha, cl, 1);
    cl_slope = pol_fit(1);
    alpha_0 = -pol_fit(2)/pol_fit(1);

end


% this function makes the flap analysis for the second section of the
% homewrk
function checkFlapEfficiency(f, p, E_ratios)
    
    flap_efficiency = zeros(1, size(E_ratios, 2));
    
    [trash, alpha_0, trash] = getAirfoilCharacteristics(f, p, -1, 0);
    
    for i = 1:size(E_ratios,2)
        efficiency = zeros(1, 40);
        for j = 1:45
            [trash, alpha(j), trash] = getAirfoilCharacteristics(f, p, (1-E_ratios(i))*10, j);
            efficiency = (alpha - alpha_0)/j;
        end
        
        efficiency
    end

end


% This function makes everything for the third section of the homework
% such a mess
function discussion()

    % values of f
    all_f = [0 0.02 0.04 0.06];
    all_p = [0.1 0.2 0.4 0.6];
    
    % two tables to hold all the obtained values
    alpha = zeros(size(all_f,2), size(all_p, 2));
    cm0 = zeros(size(all_f,2), size(all_p, 2));
    
    for f = 1:size(all_f,2)
        for p = 1:size(all_p, 2)
            [trash, alpha(f,p), cm0(f,p)] = getAirfoilCharacteristics(all_f(f)*100, all_p(p)*10, -1, 0);
        end
    end
    
    subplot(1, 2, 1);
    % jesus christ thats horrible code
    plot(all_p, alpha(1,:), all_p, alpha(2,:), all_p, alpha(3,:), all_p, alpha(4,:));
    title('zero angle of attack as a function of f and p');
    xlabel('value of p (per unit of chord)');
    ylabel('alpha 0 (º)');
    lgd = legend('f = 0.00', 'f = 0.02', 'f = 0.04', 'f = 0.06');
    lgd.Location = 'southwest';
    grid on;
    
    subplot(1, 2, 2);
    plot(all_p, cm0(1,:), all_p, cm0(2,:), all_p, cm0(3,:), all_p, cm0(4,:));
    ylim([-0.25 0.05]);
    title('Cm as a function of f and p');
    xlabel('value of p (per unit of chord)');
    ylabel('Cm');
    lgd = legend('f = 0.00', 'f = 0.02', 'f = 0.04', 'f = 0.06');
    lgd.Location = 'southwest';
    grid on;

end