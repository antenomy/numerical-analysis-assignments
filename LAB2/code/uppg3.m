load dollarkurs.mat
X = USDSEK;
N = length(X);
tt=(1:N)';



%% 3a Linjär modell


function [Q, R] = gramSchmidt(A)
    [m, n] = size(A);
    Q = zeros(m,m);
    R = zeros(m,n);
    for j=1:1:n
        y = A(:,j);
        for i=1:1:j-1
            R(i,j) = dot(Q(:,i),A(:,j));
            y = y - R(i,j)*Q(:,i);
        end
        R(j,j) = norm(y);
        Q(:,j) = y/R(j,j);
    end
end

%[testQ , testR] = gramSchmidt(A_linear);

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
    %[m, n] = size(A);
    %R_hat = R(1:n,1:n);
    %d = (Q') * b;
    %d_hat = d(1:n);
    %x = R_hat \ d_hat;
    
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
    % Er kod här...
end

function MSE = MSE(X, Y)
    N = length(X);
    MSE = (1/N)*sum((X-Y).^2);
end

b = X(:);

A_linear = [ones(N, 1), tt];

L = 500;
A_periodic = zeros(N, 4);
for t = 1:N
    A_periodic(t, :) = [1, t, sin(2*pi*t/L), cos(2*pi*t/L)];
end



%[linear_coeffs, cond_num] = leastSquares(A_linear, b);
%[linear_coeffs_QR, cond_num_QR] = leastSquaresQR(A_linear, b);

[periodic_coeffs, periodic_cond_num] = leastSquaresQR(A_periodic, b);

%disp(linear_coeffs);
%disp(linear_coeffs_QR);
%disp(condition_num);
%disp(condition_num_QR);

%linear_model = linear_coeffs(1) + linear_coeffs(2) * tt;
periodic_model = periodic_coeffs(1)  ...
    + periodic_coeffs(2) * tt ...
    + periodic_coeffs(3) * sin(2*pi*tt/L) ...
    + periodic_coeffs(4) * cos(2*pi*tt/L);
%plotSeries(X, tt, linear_model);
plotSeries(X, tt, periodic_model);
%fprintf('Linear model:\ty = %f + %f*t\n', linear_coeffs);
%fprintf('Mean squared error: %f\n', MSE(X, linear_model))


%% 3b Normalekvationerna och konditionering

%% 3c Linjär + periodisk modell

% Er kod här...


%% 3d Icke-linjär modell

% Er kod här...

