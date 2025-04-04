%% Interpolation



h = 0.25; % Får ej ändras i koden nedan
real_h = 10^-5;

%% Linjär interpolation

function [coefficient_array] = linear_interpolation(x_array, y_array)
    m = (y_arrray(2) - y_arrray(1)) / (x_array(2) - x_array(1));
    k = y_array(1) - m * x_array(1);
    coefficient_array = [m, k];
end


%% Kvadratisk interpolation

function [coefficient_array] = quadratic_interpolation(x_array, y_array)

    n_matrix = zeros(3, 3);

    for iteration = 1:length(y_array)
        n_matrix(iteration, 1) = y_array(iteration);
    end

    for i = 2:length(y_array)
        for j = 1:length(y_array)+1-i
            n_matrix(j, i) = (n_matrix(j+1, i-1)-n_matrix(j, i-1))/(x_array(j+i-1)-x_array(j)); 
        end
    end

    a = n_matrix(1, 3);
    b = n_matrix(1, 3) * (-x_array(1)-x_array(2)) + n_matrix(1, 2);
    c = n_matrix(1, 3) * (x_array(1)*x_array(2))  - n_matrix(1, 2) * x_array(1) + n_matrix(1, 1);

    coefficient_array = [a, b, c];
end

function [px,py] = piecewise_quadratic()
    for i = 1:length(x)-2
        piece_coeff = quadratic_interpolation(x(i:i+2), y(i:i+2));

    end
end

function [y] = homemade_polyval(poly_coeffs, x)
    deg = length(poly_coeffs) - 1;
    n = length(x);
    y = zeros(n, 1);
    for i=1:n
        terms = zeros(deg, 1);
        for j=deg:-1:1
            terms(j) = poly_coeffs(j) * x(i)^j;
        end
        y(i) = sum(terms);
    end
end

[t,x,y,vx,vy] = kastbana(h);
[real_t,real_x,real_y,real_vx,real_vy] = kastbana(real_h);

%coeff_linear = linear_interpolation()

coeff_quadratic = quadratic_interpolation(x(10:12), y(10:12));

plot_x = 0:h:100;
plot_y = polyval(coeff_quadratic, plot_x);
test_plot_y = homemade_polyval;

plot(x,y, "b", plot_x, plot_y, "--r");


% Calculation of max_x, max_y, down_x Quadratic 

for i = 1:length(x)-2
    if y(i) < y(i+1) && y(i+1) > y(i+2)
        max_coeff_quad = quadratic_interpolation(x(i:i+2), y(i:i+2));
        max_x_quad = -max_coeff(2)/max_coeff(1);
        max_y_quad = polyval(max_coeff, max_x);
        disp(max_x);
        disp(max_y);
        break
    end

    if y(i) == 0
        down_x_quad = 
    end
end

% Calculation of max_x, max_y, down_x Linear 

for i = 1:length(x)-2

    % Finding max x, y
    if y(i) <= y(i+1) && y(i+1) >= y(i+2)
        % max_coeff_lin = linear_interpolation(x(i:i+1), y(i:i+1));

        max_x_lin = x(i+1);
        max_x_error_lin = 0; % Need to find the two points around max_x_lin and evaluate the line at that x

        max_y_lin = y(i+1);
        max_y_error_lin = ;
    end
end

for i = 1:length(x)-1
    % Finding down x
    if y(i) == 0 
        down_x_lin = x(i);
    elseif y(i + 1) == 0 
        down_x_lin = x(i+1);
    elseif (y(i) < 0 && y(i+1) > 0) || (y(i) > 0 && y(i+1) < 0)
        down_coeff_lin = linear_interpolation(x(i:i+1), y(i:i+1));
        down_x_lin = x(i);
    end
end

disp("Linear Interpolation")
disp("Max x: ", max_x_lin, "  Max x error:", max_x_error_lin);
disp("Max y: ", max_y_lin, "  Max y error:", max_y_error_lin);
disp("Down x: ", down_x_lin, "  Down x error:", down_x_error_lin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,x,y,vx,vy]=kastbana(h)

%KASTBANA(H) beräknar banan för ett kast med en liten boll.
%
%   Dynamiken ges av en ODE som inkluderar effekten av luftmotståndet,
%
%      r'' = -g*ez-sigma*r'*|r'|/m.
%
%   Funktionen beräknar bollens position och hastighet vid
%   tidpunkter separerade med en given steglängd. Bollen kastas från
%   (X,Y)=(0,2) med hastigheten 30 m/s i 45 graders vinkel uppåt.
%
%   Syntax:
%
%   [T,X,Y,VX,VY] = KASTBANA(H)
%
%   H       - Steglängden mellan tidpunkterna.
%   T       - Vektor med tidpunkter där bollens position och hastighet beräknats.
%   X, Y    - Vektorer med bollens x- och y-koordinater vid tidpunkterna.
%   VX, VY  - Vektorer med bollens hastigheter i x- och y-led vid tidpunkterna.

%% Tennisboll, specifikationer

m = 56e-3;     % Massan (kg) = 56 gram
ra = 6.6e-2/2; % 6.6 cm in diameter

g=9.81;      % Tyngdaccelerationen (m/s^2)

rho=1.2;     % Luftens densitet (kg/m^3)
A=ra^2*pi;   % Kroppens tvärsnittsarea (m^2)
Cd=0.47;     % Luftmotståndskoefficient,
% "drag coefficient" (dimensionslös)
% Läs mer på http://en.wikipedia.org/wiki/Drag_coefficient

sigma = rho*A*Cd/2; % Totala luftmotståndet

T  = 5;      % Sluttid
v0 = 32;     % Utkasthastighet
al = pi/4;   % Utkastvinkel

% Begynnelsevärden

r0 = [0 2]';                   % Position
r1 = [v0*cos(al) v0*sin(al)]'; % Hastighet

% ODEns högerled

f = @(u) [u(3:4); -u(3:4)*norm(u(3:4),2)*sigma/m - [0;g]];  % RHS

u = [r0;r1];
U = u';
t = 0:h:T;

% Runge-Kutta 4

for tn=t(1:end-1)
    s1 = f(u);
    s2 = f(u + h/2*s1);
    s3 = f(u + h/2*s2);
    s4 = f(u + h*s3);
    u = u + h/6*(s1 + 2*s2 + 2*s3 + s4);
    U = [U; u'];
end

x  = U(:,1);
y  = U(:,2);
vx = U(:,3);
vy = U(:,4);

end

