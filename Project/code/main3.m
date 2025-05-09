format long

%omega = 19;

%S0 = @(r) cos(24*sqrt(r)).*exp(-900*r);

%aa = 1;
%xs  = 0.45;
%ys  = 0.2;
%S = @(x,y) aa*S0((x-xs).^2+(y-ys).^2);

%etaIntegrand = @(x,y, omega) S0(x.^2+y.^2).*cos(omega*x);
%eta = trapInt2D(200, etaIntegrand, omega);
%[B,Sol] = hhsolver(omega,S,400); 


%noiselevel = 0;
%gnoise = B.un+max(abs(B.un)).*randn(size(B.un)).*noiselevel;

%MAX_ALPHA = 100;
%alphas = linspace(0, 2*pi, MAX_ALPHA);
%Ic_alphas = zeros(length(alphas), 1);
%for i=1:length(alphas)
%    Ic_alphas(i) = getIc_alpha(B, alphas(i), omega, v_c, gnoise);
%end

%startGuess = getStartGuess(B, eta, omega, v_c, v_s, gnoise);

%disp(startGuess)

%nonlinear_coeffs = gaussNewton(10, length(alphas), startGuess, eta, omega, alphas, Ic_alphas);
%disp(nonlinear_coeffs)
%plotAlphas(alphas, Ic_alphas, aa, eta, v_c);
%plotFields(B, Sol, S);


%S2 = @(x,y) nonlinear_coeffs(1)*S0((x-nonlinear_coeffs(2)).^2+(y-nonlinear_coeffs(3)).^2);
%[B2,Sol2] = hhsolver(omega,S2,100); 

%plotFields(B2, Sol2, S2);

params.HHSOLVER_STEPS = 400; % How precise should hhsolver solution be
params.ALPHA_NUM = 100; % How many different angles should be tested
params.ETA_TRAPEZOID_STEPS = 800; % How many steps used in integrating for eta (for x and y domain individually)
params.GAUSS_NEWTON_RUNS = 12; % How many iterations Gauss-Newton algorithm should run
params.DERIV_STEP = 0.000001; % Step size in calculating derivatives (for I_c, I_s)

MODE = "CUSTOM";
SOURCE_NUM = 1;

OMEGA = 19;
AMPLITUDE = 1;
X_SOURCE  = 0.45;
Y_SOURCE  = 0.2;
NOISE_LEVEL = 10; % Add noise to sound boundary data


filenames = ["source1.mat", "source2.mat", "source3.mat", "source4.mat", "source5.mat"];


S0 = @(r) cos(24*sqrt(r)).*exp(-900*r);
v_c = @(x, y, alpha, omega) cos(omega.*(x.*cos(alpha) + y.*sin(alpha)));
v_s = @(x, y, alpha, omega) sin(omega.*(x*cos(alpha) + y*sin(alpha)));

switch MODE
    case "FILE" % Solve a problem from one of the files
        file = load(filenames(SOURCE_NUM));
        OMEGA = file.omega;
        disp(OMEGA)
        B = file.B;
        solveSoundProblem(B, OMEGA, S0, v_c, v_s, params);
    case "CUSTOM" % Solve custom problem
        S = @(x,y) AMPLITUDE*S0((x-X_SOURCE).^2+(y-Y_SOURCE).^2);
        [B,Sol] = hhsolver(OMEGA,S,params.HHSOLVER_STEPS);
        
        B.un = B.un+max(abs(B.un)).*randn(size(B.un)).*NOISE_LEVEL;
        [amplitude, x0, y0, eta, alphas, Ic_alphas] = solveSoundProblem(B, OMEGA, S0, v_c, v_s, params);
        plotAlphas(alphas, Ic_alphas, AMPLITUDE, X_SOURCE, Y_SOURCE, OMEGA, eta, v_c)
    case "NOISE"
        S = @(x,y) AMPLITUDE*S0((x-X_SOURCE).^2+(y-Y_SOURCE).^2);
        [B,Sol] = hhsolver(OMEGA,S,params.HHSOLVER_STEPS);
        levels = 0; %[1e-2, 1e-1, 0.5, 1.0];
        figure;
        for k=1:length(levels)
            lvl = levels(k);
            gnoise = B.un + max(abs(B.un)) * randn(size(B.un)) * lvl;
            %subplot(2,2,k);
            plot(B.s, B.un,    'b-'); hold on;
            %plot(B.s, gnoise,'r.');
            hold off;
            title(['g = ',num2str(lvl)]);
            xlabel('s'); ylabel('g(s)');
            axis tight; grid on;
        end
end




function [amplitude, x0, y0, eta, alphas, Ic_alphas] = solveSoundProblem(BB, omega, S0Func, v_c, v_s, params)
    etaIntegrand = @(x,y, omega) S0Func(x.^2+y.^2).*cos(omega*x);
    eta = trapInt2D(params.ETA_TRAPEZOID_STEPS, etaIntegrand, omega);
    disp(eta)
    
    alphas = linspace(0, 2*pi, params.ALPHA_NUM);
    Ic_alphas = zeros(length(alphas), 1);
    for i=1:length(alphas)
        Ic_alphas(i) = getIc_alpha(BB, alphas(i), omega, v_c);
    end

    startGuess = getStartGuess(BB, eta, omega, v_c, v_s, params.DERIV_STEP);
    problem_coeffs = gaussNewton(params.GAUSS_NEWTON_RUNS, length(alphas), startGuess, eta, omega, alphas, Ic_alphas);
    
    amplitude = problem_coeffs(1);
    x0 = problem_coeffs(2);
    y0 = problem_coeffs(3);
    
    disp(startGuess)

    disp(['Omega: ', num2str(omega)])
    disp(['Amplitude: ', num2str(amplitude, 15)])
    disp(['x0: ', num2str(x0, 15)])
    disp(['y0: ', num2str(y0, 15)])   

    SFunc = @(x,y) amplitude*S0Func((x-x0).^2+(y-y0).^2);
    [estimatedB, estimatedSol] = hhsolver(omega,SFunc,params.HHSOLVER_STEPS); 
    plotFields(BB, estimatedB, estimatedSol, SFunc);

    
end

function I = trapInt2D(n, func, omega)
    a = -4;
    b = 4;
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

function x = gaussNewton(k, N, startGuess, eta, omega, alphas, Ic_alphas)
    %x = [1.05, 0.4, 0.25]';
    x = startGuess';
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

function deriv = threePointDeriv(func, alpha, h)
    % No need to wrap around alpha values inside [0, 2pi] since
    %   all associated functions are already periodic such that
    %   we only care about the equivalence class.
    f_left = func(alpha-h);
    f_right = func(alpha+h);

    deriv = (f_right - f_left) ./ (2*h);
end

function Ic_alpha = getIc_alpha(B, alpha, omega, vc_func)
    % Calculate v_c test wave vals for boundary points
    
    vc_vals = vc_func(B.x, B.y, alpha, omega);

    % Integrand for I_c(alpha)
    integrand = B.un .* vc_vals;
    Ic_alpha = numInt(B.s, integrand);
end

function Is_alpha = getIs_alpha(B, alpha, omega, vs_func)
    % Calculate v_s test wave vals for boundary points
    
    
    vs_vals = vs_func(B.x, B.y, alpha, omega);

    % Integrand for I_s(alpha)
    integrand = B.un .* vs_vals;
    Is_alpha = numInt(B.s, integrand);
end

function startGuess = getStartGuess(B, eta, omega, vc_func, vs_func, h)
    IcHandle = @(a) getIc_alpha(B, a, omega, vc_func);
    IsHandle = @(a) getIs_alpha(B, a, omega, vs_func);

    % Ic values
    Ic_0 = IcHandle(0);
    dIc_0 = threePointDeriv(IcHandle, 0, h);
    dIc_pi2 = threePointDeriv(IcHandle, pi/2, h);

    % Is values
    Is_0 = IsHandle(0);
    Is_pi2 = IsHandle(pi/2);

    a_0 = sqrt(Ic_0.^2 + Is_0.^2) / eta;
    x_0 = dIc_pi2 / (omega * Is_pi2);
    y_0 = -1 * dIc_0 / (omega * Is_0);
    
    startGuess = [a_0, x_0, y_0];
end

function plotFields(BB, estimatedB, Sol, sourceModel)
    figure;

    surf(Sol.x,Sol.y,Sol.u)
    view(0,90)
    shading interp
    axis equal

    figure;

    surf(Sol.x,Sol.y,Sol.u)
    shading interp
    light                    
    lighting gouraud       % Tar lång tid när N är stor
    view(-60,60)
    
    figure;
    %surf(Sol.x,Sol.y,Sol.S)
    surf(Sol.x,Sol.y,sourceModel(Sol.x,Sol.y))
    shading interp
    light
    lighting gouraud       % Tar lång tid när N är stor
    colormap('autumn')
    material shiny
    view(-60,60)
    
    figure;
    contour(Sol.x,Sol.y,Sol.u,20)
    axis equal
    hold on
    plot(BB.x,BB.y,'k-','LineWidth',2)
    [c,hnd]=contour(Sol.x,Sol.y,sourceModel(Sol.x,Sol.y),10); %Sol.S,10);
    set(hnd,'Color','k','LineWidth',1.5)
    hold off
    axis off

    figure;
    plot3(BB.x,BB.y,BB.un,'b-','LineWidth',2) % Actual data
    hold on
    plot3(estimatedB.x,estimatedB.y,estimatedB.un,'r-','LineWidth',2) % Estimated data
    contour(Sol.x,Sol.y,Sol.u,20)
    plot3(BB.x,BB.y,zeros(size(BB.x)),'k-','LineWidth',2)
    hold off
        
    figure;
    mesh(Sol.x,Sol.y,Sol.u)
end

function plotAlphas(alphas, Ic_alphas, aa, xs, ys, omega, eta, v_c)
    % Compute the right-hand side: aa * eta * v_c
    vc_vals = zeros(length(alphas), 1);
    for i = 1:length(alphas)
        vc_vals(i) = v_c(xs, ys, alphas(i), omega);
    end
    rhs = aa * eta * vc_vals;
    figure;
    plot(alphas, Ic_alphas, 'b-'); 
    hold on;
    grid on;
    plot(alphas, rhs, 'r--');
    hold off;
end

