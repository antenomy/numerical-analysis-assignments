%% 1b -- Visualisering


file = 'eiffel1.mat';
load(file)

matrix_height = height(A);
b = zeros(matrix_height, 1);

%[V, D] = eig(A);
%lambdas = diag(D);

%[sorted_lambdas, index] = sort(lambdas, "descend");
%sorted_V = V(:, index);

%disp('Sorted Eigenvalues:');
%disp(sorted_lambdas);

%disp('Sorted Eigenvectors:');
%disp(sorted_V);

%frequencies = sqrt(sorted_lambdas) / (2 * pi);
%mode_index = size(sorted_V,2); 

for i = 1:1
    break % Remove when finished task
    mode = sorted_V(:, mode_index-(i-1));

    xEigen = xnod + 5 * mode(1:2:end);
    yEigen = ynod + 5 * mode(2:2:end);

    figure(i) 
    trussanim(xnod, ynod, bars, mode);
    hold on
    title(['Eigenmode ' num2str(i)])
    fprintf('Frequency:%s\n', frequencies(mode_index-(i-1)));
    hold off
end



%% 1c -- Beräkning av största och minsta egenvärdena

file = {'eiffel1.mat', 'eiffel2.mat', 'eiffel3.mat', 'eiffel4.mat'};

for i = 1:length(file)
    fprintf('\neiffel%.0f\n', i);
    tau = 10e-10;
    [mu, iter] = potens(A,tau);
    fprintf('[Power Iteration]\n');
    fprintf('Eigenvalue: %f     Iterations: %.0f\n', mu, iter);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    [mu, iter] = inverspotens(A,tau);
    fprintf('[Inverse Power Iteration]\n');
    fprintf('Eigenvalue: %f     Iterations: %.0f\n', mu, iter);
end



%% 1d -- Ber�kning av andra egenv�rden

% Er kod h�r...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [mu, iter] = potens(A,tau)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Indata:
    %
    %  A  - matrisen (kvadratisk)
    %  tau - feltolerans (sk�l�r)
    %
    %  Utdata:
    %
    %  mu - st�rsta egenv�rdet till A (skal�r)
    %  iter - antal iterationer som anv�nts (skal�r)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = ones(height(A), 1);
    iter = 0;
    while (0==0)
        iter = iter + 1;

        u = x/norm(x);      % normalise vector
        x = A*u;           % power step
        mu = u'*x;         % Rayleigh Quotient

        diff = norm((A - mu*eye(size(A)))*u);
        if (diff < tau)
            break
        end
        if iter >= 1000
            disp("Timeout")
            break
        end
    end
    u = x/norm(x);
end

function [mu, iter] = inverspotens(A,tau)

    %  Indata:
    %
    %  A  - matrisen (kvadratisk)
    %  tau - feltolerans (sk�l�r)
    %
    %  Utdata:
    %
    %  mu - minsta egenv�rdet till A (skal�r)
    %  iter - antal iterationer som anv�nts (skal�r)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = ones(height(A), 1);
    iter = 0;
    while (0==0)
        iter = iter + 1;
        
        u = x/norm(x);      % normalise vector
        x = A\u;           % power step
        mu = u'*x;         % Rayleigh Quotient
        

        diff = norm((A - (1/mu)*eye(size(A)))*u);

        if (diff < tau)
            break
        end
        if iter >= 100
            disp("Timeout")
            break
        end
    end
    mu = 1/mu; 
    u = x/norm(x);
end