omega = 19;

S0 = @(r) cos(24*sqrt(r)).*exp(-900*r);
v_c = @(x, y, alpha, omega) cos(omega*(x*cos(alpha) + y*sin(alpha)));

aa = 1;
xs  = 0.45;
ys  = 0.2;
S = @(x,y) aa*S0((x-xs).^2+(y-ys).^2);

nuIntegrand = @(x,y, omega) S0(x.^2+y.^2).*cos(omega*x);

function I = trapInt2D(n, func, omega)
    a = -2;
    b = 2;
    h = (b - a) / n;
    
    x = linspace(a, b, n+1);
    y = linspace(a, b, n+1);

    I = 0;
    for j1=1:n+1
        for j2 = 1:n+1
            % Trapezoid rule weight logic
            w = 1;
            if j1 == 1 || j1 == n+1
                w = w / 2;
            end
            if j2 == 1 || j2 == n+1
                w = w / 2;
            end
            
            % function eval
            I = I + w*func(x(j1), y(j2), omega);
        end
    end

    I = I * h^2;
end

function I = numInt(s, integrand)
    n = length(s);
    I = 0;
    for i=1:n-1
        
        h = s(i+1) - s(i);

        I = I + (h/2).*(integrand(i)+integrand(i+1));
    end
end

eta = trapInt2D(100, nuIntegrand, omega);
disp(eta);

[B,Sol] = hhsolver(omega,S,100); 



function r = func_r(x, N, eta, omega, alphas)
    r = zeros(N, 1);
    for j = 1:N
        r(t) = x(1)*eta*cos(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j)));
    end
end

function Dr = derivative_r(x, N, eta, omega, alphas)
    Dr = zeros(N, 3);
    for j = 1:N
        Dr(t, :) = [eta*cos(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j))), ...
            -x(1)*omega*cos(alphas(j))*eta*sin(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j))), ...
            -x(1)*omega*sin(alphas(j))*eta*sin(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j)))];
        
        [1, t, sin(2*pi*t/x(5)), cos(2*pi*t/x(5)), ...
        (2*pi*t/(x(5)^2))*(x(4)*sin(2*pi*t/x(5))-x(3)*cos(2*pi*t/x(5)))];
    end
end

function x = gaussNewton(k, N, eta, omega, alphas)
    x = [0.9, 0.5, 0.3]';
    for i=1:k
        A = derivative_r(x, N, eta, omega, alphas);
        v = (A')*A \ -(A')*func_r(x, N, eta, omega, alphas);
        x = x + v;
    end
end

disp(B)
MAX_ALPHA = 10;
for i=1:MAX_ALPHA
    alpha = (i/MAX_ALPHA).*2*pi;

    % Calculate v_c test wave vals for boundary points
    vc_vals = v_c(B.x, B.y, alpha, omega);

    % Integrand for I_c(alpha)
    integrand = B.un .* vc_vals;
    I_c_alpha = numInt(B.s, integrand);

    disp(I_c_alpha)
    % To-do finish 3b (GN method to estimate values for a, x_0 and y_0)
    nonlinear_coeffs = gaussNewton();

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
