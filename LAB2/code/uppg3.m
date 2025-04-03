load dollarkurs.mat
X = USDSEK;
N = length(X);
tt=(1:N)';
b = X(:);



%% 3a Linjär modell

A_linear = [ones(N, 1), tt];
[linear_coeffs_QR, cond_num_QR] = leastSquaresQR(A_linear, b);
linear_model = linear_coeffs_QR(1) + linear_coeffs_QR(2) * tt;

fprintf('Linear model:\ty = %f + %f*t\n', linear_coeffs_QR);
fprintf('Mean squared error: %f\n', MSE(X, linear_model))

plotSeries(X, tt, linear_model);


%% 3b Normalekvationerna och konditionering

A_linear = [ones(N, 1), tt];
[linear_coeffs, cond_num] = leastSquares(A_linear, b);
[linear_coeffs_QR, cond_num_QR] = leastSquaresQR(A_linear, b);
disp(linear_coeffs);
disp(linear_coeffs_QR);
disp(cond_num);
disp(cond_num_QR);


%% 3c Linjär + periodisk modell

L = 500;
A_periodic = zeros(N, 4);
for t = 1:N
    A_periodic(t, :) = [1, t, sin(2*pi*t/L), cos(2*pi*t/L)];
end

[periodic_coeffs, periodic_cond_num] = leastSquaresQR(A_periodic, b);

periodic_model = periodic_coeffs(1)  ...
    + periodic_coeffs(2) * tt ...
    + periodic_coeffs(3) * sin(2*pi*tt/L) ...
    + periodic_coeffs(4) * cos(2*pi*tt/L);

fprintf('Periodic model:\ty = %f + %f*t + %f*sin(2*pi*t/L) + %f*cos(2*pi*t/L), L = 500\n', periodic_coeffs);
fprintf('Mean squared error: %f\n', MSE(X, periodic_model))

plotSeries(X, tt, periodic_model);



%% 3d Icke-linjär modell



nonlin_coeffs = gaussNewton(10, N, b);
nonlin_model = nonlin_coeffs(1)  ...
    + nonlin_coeffs(2) * tt ...
    + nonlin_coeffs(3) * sin(2*pi*tt/nonlin_coeffs(5)) ...
    + nonlin_coeffs(4) * cos(2*pi*tt/nonlin_coeffs(5));

fprintf('Nonlinear model:\ty = %f + %f*t + %f*sin(2*pi*t/L) + %f*cos(2*pi*t/L), L = %f\n', nonlin_coeffs);
fprintf('Mean squared error: %f\n', MSE(X, nonlin_model))

plotSeries(X, tt, nonlin_model);

% FUNCTIONS

function [x, C] = leastSquares(A, b)
    ATA = (A')*A;
    x = ATA \ (A')*b;
    C = cond(ATA);
end

function [x, C] = leastSquaresQR(A, b)
    [Q, R] = qr(A, 0); % Second parameter indicates we want the reduced factorization, 
    %   which saves computing power since we onlt need Q to have the basis elements
    %   for the basis for the column space of A (n=2) rather than the entire space 
    %   which is isomorphic to R^(m*m) (m=big number).

    x = R \ (Q') * b; %  Back substitution for QR variation of normal equations
    C = cond(R);
end

function plotSeries(X, tt, linear_model)
    ts = timeseries(X);
    ts.Name = "USD-SEK";
    ts.TimeInfo.Units = "days";
    ts.TimeInfo.StartDate = '2009-01-01';
    ts.TimeInfo.Format = 'mmmm dd, yyyy';
    plot(ts);
    hold on;

    plot(tt, linear_model, 'r-', 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('USD-SEK');
    hold off;
    waitfor(gcf);
end

function MSE = MSE(X, Y)
    N = length(X);
    MSE = (1/N)*sum((X-Y).^2);
end

function r = func_r(x, N, b)
    r = zeros(N, 1);
    for t = 1:N
        r(t) = x(1) + x(2)*t + x(3)*sin(2*pi*t/x(5)) + x(4)*cos(2*pi*t/x(5)) - b(t);
    end
end

function Dr = derivative_r(x, N)
    Dr = zeros(N, 5);
    for t = 1:N
        Dr(t, :) = [1, t, sin(2*pi*t/x(5)), cos(2*pi*t/x(5)), ...
        (2*pi*t/(x(5)^2))*(x(4)*sin(2*pi*t/x(5))-x(3)*cos(2*pi*t/x(5)))];
    end
end

function x = gaussNewton(k, N, b)
    x = [8, -0.002, 0.3, 0.4, 500]';
    for i=1:k
        A = derivative_r(x, N);
        v = (A')*A \ -(A')*func_r(x, N, b);
        x = x + v;
    end
end