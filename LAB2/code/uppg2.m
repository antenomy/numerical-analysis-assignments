%% Interpolation



h = 0.25; % Får ej ändras i koden nedan
real_h = 10e-5;

%% Linjär interpolation

function [coefficient_array] = linear_interpolation(x_array, y_array)
    m = (y_array(2) - y_array(1)) / (x_array(2) - x_array(1));
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


function [y] = homemade_polyval(poly_coeffs, x)
    deg = length(poly_coeffs) - 1;
    n = length(x);
    y = zeros(1, n);
    for i=1:n
        terms = zeros(deg+1, 1);
        for j=1:deg+1
            terms(j) = poly_coeffs(j) * x(i)^(deg-j+1);
        end
        y(1,i) = sum(terms);
    end
end

[t,x,y,vx,vy] = kastbana(h);
[real_t,real_x,real_y,real_vx,real_vy] = kastbana(real_h);

disp(size(x))
disp(size(y))
plot(x, y, 'bo-', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b')
hold on;
disp(size(real_x))
disp(size(real_y))
plot(real_x, real_y, '--g')
hold on;


% Calculation of max_x, max_y, down_x Quadratic 

function x_max_help = calculate_max_x_quad_helper(coeff)
    x_max_help = -coeff(2)/(2*coeff(1));
end

function [x_max, y_max] = calculate_max_x_quad(x, y)
    for i = 1:length(x)-2
        if y(i) <= y(i+1) && y(i+1) >= y(i+2) % Finds local maximum
            
    
            if mod(i, 2) == 1 % One interpolating polynomial to check
                max_coeff_quad = quadratic_interpolation(x(i-1:i+1), y(i-1:i+1));
                x_max = calculate_max_x_quad_helper(max_coeff_quad);
                y_max = homemade_polyval(max_coeff_quad, x_max);

                plot_x = 0:0.2:100;
                plot_y = homemade_polyval(max_coeff_quad, plot_x);
                plot(plot_x, plot_y, "--r");
                hold on;
                grid on;
                return;
            else % Two interpolating polynomials to check
                max_coeff_quad_1 = quadratic_interpolation(x(i-2:i), y(i-2:i));
                max_coeff_quad_2 = quadratic_interpolation(x(i:i+2), y(i:i+2));

                plot_x = 0:0.2:100;
                plot_y = homemade_polyval(max_coeff_quad_1, plot_x);
                plot(plot_x, plot_y, "--r");
                hold on;
                plot_y = homemade_polyval(max_coeff_quad_2, plot_x);
                
                plot(plot_x, plot_y, "--r");
                grid on;

                x_max_1 = calculate_max_x_quad_helper(max_coeff_quad_1);
                x_max_2 = calculate_max_x_quad_helper(max_coeff_quad_2);

                in_boundary_1 = (x_max_1 <= x(i-2) && x_max_1 >= x(i)) || (x_max_1 >= x(i-2) && x_max_1 <= x(i));
                in_boundary_2 = (x_max_2 <= x(i) && x_max_2 >= x(i+2)) || (x_max_2 >= x(i) && x_max_2 <= x(i+2));
                
                if in_boundary_1 && in_boundary_2
                    x_max = max(x_max_1, x_max_2);
                elseif in_boundary_1
                    x_max = x_max_1;
                elseif in_boundary_2
                    x_max = x_max_2;
                else
                    x_max = x(i);
                end
                y_max = y(i);
                return;
            end
        end
    end
    x_max = 0;
    y_max = 0;
    return;
end


function down_x_quad = calculate_down_x_quad(x, y)
    for i = 1:2:length(x)-2
        % Finding down x
        if y(i) == 0 
            down_x_quad = x(i);
            return;
        elseif y(i + 1) == 0 
            down_x_quad = x(i+1);
            return;
        elseif y(i + 2) == 0 
            down_x_quad = x(i+2);
            return;
        elseif y(i)*y(i+2) < 0
            down_coeff_quad = quadratic_interpolation(x(i:i+2), y(i:i+2));
            
            down_x_quad_array = quadratic_formula(down_coeff_quad);

            if (down_x_quad_array(1) <= x(i) && down_x_quad_array(1) >= x(i+2)) || (down_x_quad_array(1) >= x(i) && down_x_quad_array(1) <= x(i+2))
                down_x_quad = down_x_quad_array(1);
            elseif (down_x_quad_array(2) <= x(i) && down_x_quad_array(2) >= x(i+2)) || (down_x_quad_array(2) >= x(i) && down_x_quad_array(2) <= x(i+2))
                down_x_quad = down_x_quad_array(2);
            return;
            end
        end
    end

    down_x_quad = 0;
    return;
end

function x_array = quadratic_formula(coeff_array)
    x_array = zeros(1, 2);
    % Extract coefficients
    a = coeff_array(1);
    b = coeff_array(2);
    c = coeff_array(3);
    
    % Compute discriminant
    discriminant = b^2 - 4*a*c;
    
    % Compute root(s)
    x_array(1) = (-b +  1 * sqrt(discriminant)) / (2*a);
    x_array(2) = (-b + -1 * sqrt(discriminant)) / (2*a);
end

[max_x_quad, max_y_quad] = calculate_max_x_quad(x,y);
[max_x_quad_real, max_y_quad_real] = calculate_max_x_quad(real_x,real_y);

max_x_quad_error = abs(max_x_quad - max_x_quad_real);
max_y_quad_error = abs(max_y_quad - max_y_quad_real);

down_x_quad = calculate_down_x_quad(x, y);
down_x_quad_real = calculate_down_x_quad(real_x, real_y);

down_x_quad_error = abs(down_x_quad - down_x_quad_real);

disp("Quadratic Interpolation")
fprintf('Max x: %f Max x error: %f \n', max_x_quad, max_x_quad_error);
fprintf('Max y: %f Max y error: %f \n\n', max_y_quad, max_y_quad_error);

%disp("Max y: ", max_y_lin, "  Max y error:", max_y_error_lin);
fprintf('Down x: %f Down x error: %f \n\n', down_x_quad, down_x_quad_error);



% Calculation of max_x, max_y, down_x Linear 

% Finding max x, y
function [x_max, y_max] = calculate_max_lin(x, y)
    n = length(x);

    i_save = 1;

    for i = 1:n

        if y(i_save) < y(i)
            i_save = i;
        end
    end

    x_max = x(i_save);
    y_max = y(i_save);
    return;
end

function down_x_lin = calculate_down_x_lin(x, y)
    for i = 1:length(x)-1
        % Finding down x
        if y(i) == 0 
            down_x_lin = x(i);
            return;
        elseif y(i + 1) == 0 
            down_x_lin = x(i+1);
            return;
        elseif y(i)*y(i+1) < 0
            down_coeff_lin = linear_interpolation(x(i:i+1), y(i:i+1));
            down_x_lin = -down_coeff_lin(2)/down_coeff_lin(1);
            return;
        end
    end

    down_x_lin = 0;
    return;
end

[max_x_lin, max_y_lin] = calculate_max_lin(x, y);
[max_x_lin_real, max_y_lin_real] = calculate_max_lin(real_x, real_y);

max_x_lin_error = max_x_lin - max_x_lin_real;
max_y_lin_error = max_y_lin - max_y_lin_real;

down_x_lin = calculate_down_x_lin(x, y);
down_x_lin_real = calculate_down_x_lin(real_x, real_y);

down_x_error_lin = abs(down_x_lin - down_x_lin_real);

disp("Linear Interpolation")
fprintf('Max x: %f Max x error: %f \n', max_x_lin, max_x_lin_error);
fprintf('Max y: %f Max y error: %f \n\n', max_y_lin, max_y_lin_error);
fprintf('Down x: %f Down x error: %f \n\n', down_x_lin, down_x_error_lin);

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

