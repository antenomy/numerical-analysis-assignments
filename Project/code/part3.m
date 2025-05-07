format long

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

function r = func_r(x, N, eta, omega, alphas, Ic_alphas)
    r = zeros(N, 1);
    for j = 1:N
        r(j) = x(1)*eta*cos(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j))) - Ic_alphas(j);
    end
end

function Dr = derivative_r(x, N, eta, omega, alphas)
    Dr = zeros(N, 3);
    for j = 1:N
        Dr(j, :) = [eta*cos(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j))), ...
            -x(1)*omega*cos(alphas(j))*eta*sin(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j))), ...
            -x(1)*omega*sin(alphas(j))*eta*sin(omega*x(2)*cos(alphas(j)) + omega*x(3)*sin(alphas(j)))];
    end
end

function x = gaussNewton(k, N, eta, omega, alphas, Ic_alphas)
    x = [1.05, 0.4, 0.25]';
    for i=1:k
        A = derivative_r(x, N, eta, omega, alphas);
        cond_num = cond(A' * A);
        if cond_num > 1e10
            warning('Matrix (A'' * A) is ill-conditioned. Results may be inaccurate.');
        end
        v = (A')*A \ -(A')*func_r(x, N, eta, omega, alphas, Ic_alphas);
        x = x + v;
    end
end

eta = trapInt2D(200, nuIntegrand, omega);
disp(eta);

[B,Sol] = hhsolver(omega,S,100); 

MAX_ALPHA = 100;
alphas = linspace(0, 2*pi, MAX_ALPHA);
Ic_alphas = zeros(length(alphas), 1);

noiselevel = 3;
gnoise = B.un+max(abs(B.un))*randn(size(B.un))*noiselevel;

for i=1:length(alphas)
    % Calculate v_c test wave vals for boundary points
    vc_vals_i = v_c(B.x, B.y, alphas(i), omega);

    % Integrand for I_c(alpha)
    integrand_i = gnoise .* vc_vals_i;
    Ic_alpha_i = numInt(B.s, integrand_i);
    Ic_alphas(i) = Ic_alpha_i;

    
end


disp(Ic_alphas)
% To-do finish 3b (GN method to estimate values for a, x_0 and y_0)
nonlinear_coeffs = gaussNewton(10, length(alphas), eta, omega, alphas, Ic_alphas);
disp(nonlinear_coeffs)


% Compute the right-hand side: aa * eta * v_c
vc_vals = zeros(length(alphas), 1);
for i = 1:length(alphas)
    vc_vals(i) = v_c(xs, ys, alphas(i), omega);
end
rhs = aa * eta * vc_vals;

figure(1);
plot(alphas, Ic_alphas, 'b-'); 
hold on;
grid on;
plot(alphas, rhs, 'r--');
hold off;

figure(2)

surf(Sol.x,Sol.y,Sol.u)
shading interp
light                    
lighting gouraud       % Tar lång tid när N är stor
view(-60,60)

figure(3)
%surf(Sol.x,Sol.y,Sol.S)
surf(Sol.x,Sol.y,S(Sol.x,Sol.y))
shading interp
light
lighting gouraud       % Tar lång tid när N är stor
colormap('autumn')
material shiny
view(-60,60)

figure(4)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off


figure(5)
mesh(Sol.x,Sol.y,Sol.u)

S2 = @(x,y) nonlinear_coeffs(1)*S0((x-nonlinear_coeffs(2)).^2+(y-nonlinear_coeffs(3)).^2);
[B2,Sol2] = hhsolver(omega,S2,100); 

figure(6)

surf(Sol2.x,Sol2.y,Sol2.u)
shading interp
light                    
lighting gouraud       % Tar lång tid när N är stor
view(-60,60)

figure(7)
%surf(Sol.x,Sol.y,Sol.S)
surf(Sol2.x,Sol2.y,S2(Sol2.x,Sol2.y))
shading interp
light
lighting gouraud       % Tar lång tid när N är stor
colormap('autumn')
material shiny
view(-60,60)

figure(8)
contour(Sol2.x,Sol2.y,Sol2.u,20)
axis equal
hold on
plot(B2.x,B2.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol2.x,Sol2.y,S2(Sol2.x,Sol2.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off


figure(9)
mesh(Sol2.x,Sol2.y,Sol2.u)