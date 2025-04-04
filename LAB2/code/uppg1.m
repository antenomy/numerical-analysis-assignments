%% 1b -- Visualisering


file = 'eiffel1.mat';
load(file)

matrix_height = height(A);
b = zeros(matrix_height, 1);

[V, D] = eig(A);
lambdas = diag(D);

[sorted_lambdas, index] = sort(lambdas, "descend");
sorted_V = V(:, index);

disp('Sorted Eigenvalues:');
disp(sorted_lambdas);

disp('Sorted Eigenvectors:');
disp(sorted_V);

frequencies = sqrt(sorted_lambdas) / (2 * pi);
mode_index = size(sorted_V,2); 

for i = 1:4
    mode = sorted_V(:, mode_index-(i-1));

    xEigen = xnod + 5 * mode(1:2:end);
    yEigen = ynod + 5 * mode(2:2:end);

    figure(i) 
    trussplot(xnod, ynod, bars, 'b');
    hold on;
    trussplot(xEigen, yEigen, bars, 'r');
    hold on;
    %trussanim(xnod, ynod, bars, mode);
    %hold on;
    title(['Eigenmode ' num2str(i)])
    fprintf('Frequency:%s\n', frequencies(mode_index-(i-1)));
    hold off
end



%% 1c -- Beräkning av största och minsta egenvärdena

file = {'eiffel1.mat', 'eiffel2.mat', 'eiffel3.mat', 'eiffel4.mat'};

T = zeros(4, 5);

for i = 1:length(file)-1
    load(file{i})
    
    %A = sparse(A);
    fprintf('\neiffel%.0f\n', i);
    tau = 10e-10;

    tic;
    [mu_pot, iter_pot] = potens(A,tau);
    potens_time = toc;
    tic;
    [mu_invpot, iter_invpot] = inverspotens(A,tau);
    inverspotens_time = toc;
    tic;
    [V, D] = eig(A);
    eig_time = toc;
    %[V, D] = eigs(A);
    lambdas = diag(D);
    
    [sorted_lambdas, index] = sort(lambdas, "descend");
    sorted_V = V(:, index);

    largest_eigenvalue = sorted_lambdas(1);
    second_largest_eigenvalue = sorted_lambdas(2);
    smallest_eigenvalue = sorted_lambdas(end);
    second_smallest_eigenvalue = sorted_lambdas(end - 1);

    ratio_largest = second_largest_eigenvalue / largest_eigenvalue;
    ratio_smallest = smallest_eigenvalue / second_smallest_eigenvalue;

    T(i, 1) = largest_eigenvalue;
    T(i, 2) = iter_pot;
    T(i, 3) = ratio_largest;
    T(i, 4) = smallest_eigenvalue;
    T(i, 5) = iter_invpot;
    T(i, 6) = ratio_smallest;

    

    fprintf('[Power Iteration]\n');
    fprintf('Eigenvalue: %f     Iterations: %.0f\n', mu_pot, iter_pot);
    fprintf('Largest Eigenvalue (eig): %f\n', largest_eigenvalue);

    
    fprintf('[Inverse Power Iteration]\n');
    fprintf('Eigenvalue: %f     Iterations: %.0f\n', mu_invpot, iter_invpot);
    fprintf('Smallest Eigenvalue (Exact): %f\n\n', smallest_eigenvalue);

    fprintf('Power Iteration time: %f\n', potens_time);
    fprintf('Inverse Power Iteration time: %f\n', inverspotens_time);
    fprintf('Eig time: %f\n', eig_time);

end

tab=array2table(T,'VariableNames',{'Största egenvärdet' '# iter potens' 'lambda2/lambda1' 'Minsta egenvärdet' '# iter inverse potens' 'lambda_n/lambda_{n-1}'},'RowNames',{'eiffel1' 'eiffel2' 'eiffel3' 'eiffel4'});
disp(tab);


%% 1d -- Ber�kning av andra egenv�rden

file = 'eiffel1.mat';
load(file)


%[mu, iter] = inverspotensmedskift(A,tau);



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
        if iter >= 10000
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
        %x = A\u;           % power step
        [L, U] = lu(A);    % LU factorisation in power step
        c = L \ u;
        x = U \ c;
        mu = u'*x;         % Rayleigh Quotient
        

        diff = norm((A - (1/mu)*eye(size(A)))*u);

        if (diff < tau)
            break
        end
        if iter >= 10000
            disp("Timeout")
            break
        end
    end
    mu = 1/mu; 
    u = x/norm(x);
end


function [mu, iter] = inverspotensmedskift(A, tau, s)

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
    x = ones(height(A), 1)./2;
    iter = 0;
    As = A-s*eye(size(A));
    while (0==0)
        iter = iter + 1;
        
        u = x/norm(x);      % normalise vector
        x = As\u;           % power step
        mu = u'*x;         % Rayleigh Quotient
        

        diff = norm((A - (1/mu + s)*eye(size(A)))*u);

        if (diff < tau)
            break
        end
        if iter >= 10000
            disp("Timeout")
            break
        end
    end
    mu = 1/mu+s; 
    u = x/norm(x);
end

disp("test");

[mu_8, iter_8] = inverspotensmedskift(A, 10e-10, 8);
[mu_55, iter_55] = inverspotensmedskift(A, 10e-10, 55);
[mu_67, iter_67] = inverspotensmedskift(A, 10e-10, 67);

fprintf('\n[Inverse Power Iteration with Shift]\n');
fprintf('Shift: 8     Eigenvalue: %.10f     Iterations: %d\n', mu_8, iter_8);
fprintf('Shift: 55    Eigenvalue: %.10f     Iterations: %d\n', mu_55, iter_55);
fprintf('Shift: 67    Eigenvalue: %.10f     Iterations: %d\n', mu_67, iter_67);