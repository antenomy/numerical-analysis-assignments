%% Assignment a)

load('eiffel1.mat');
figure(1)

matrix_size = height(A);
b = zeros(matrix_size,1); 
bel_index = round(matrix_size/6)*2 - 1;
b(bel_index) = 1;

x = A \ b;

xbel = xnod + x(1:2:end); 
ybel = ynod + x(2:2:end);

trussplot(xnod,ynod,bars, 'b')
hold on
trussplot(xbel,ybel,bars, 'r'); 
hold on
plot(xnod(bel_index), ynod(bel_index), 'y*', 'MarkerSize', 15);
hold off





%% Assignment b)

figure(2)

unkown_variables = zeros(1, 4);
compute_times = zeros(1, 4);
file = {'eiffel1.mat', 'eiffel2.mat', 'eiffel3.mat', 'eiffel4.mat'};

for iter = 1:4
    load(file{iter})
    matrix_size = height(A);
    unkown_variables(iter) = matrix_size;
    b = randn(matrix_size,1);
    total_time = 0;
    
    for iter2 = 1:20

        tic;
        x = A \ b;
        total_time = total_time + toc;

    end
    
    compute_times(iter) = total_time / 20;
end

loglog(unkown_variables, compute_times, '-o')

xlabel('Unkown Variables');
ylabel('Compute Times (s)');
grid on





%% Assignment c)

load('eiffel1.mat')
figure(3)

[min_index, max_index] = kanslighet(A, 1);

disp(['Min index: ', num2str(min_index)]);
disp(['Max index: ', num2str(max_index)]);

trussplot(xnod,ynod,bars, 'b'); 
hold on
plot(xnod(max_index), ynod(max_index), 'y*', 'MarkerSize', 15);
plot(xnod(min_index), ynod(min_index), 'yo', 'MarkerSize', 15);
hold off





%% Assignment d) - Time Table

tab=array2table(T,'VariableNames',{'Naiv' 'LU' 'Gles' 'Gles LU'},'RowNames',{'eiffel1' 'eiffel2' 'eiffel3' 'eiffel4'});
disp(tab);

function [jmin,jmax]=kanslighet(A,metod)

    matrix_size = length(A);

    Tj_max = 0;
    Tj_min = 0;

    for iter = 1:(matrix_size/2)
        b = zeros(matrix_size,1); 
        b(iter) = -1;

        if metod == 1
            x = A \ b;
        elseif metod == 2
            [L, U] = lu(A);
            c = L \ b;
            x = U \ c;
        end

        Tj = norm(x);

        if iter == 1
            Tj_max = Tj;
            jmax = 1;
        elseif Tj > Tj_max
            jmax = iter;
        end
        
        if iter == 1
            Tj_min = Tj;
            jmin = 1;
        elseif Tj < Tj_min
            jmin = iter;
        end        
    end 
end

T = zeros(4, 4);

for iter = 1:2
    load(file{iter})
    
    tic;
    kanslighet(A, 1);
    T(iter, 1) = toc
    
    tic;
    kanslighet(A, 2);
    T(iter, 2) = toc

    A = sparse(A);

    tic;
    kanslighet(A, 1);
    T(iter, 3) = toc
    
    tic;
    kanslighet(A, 2);
    T(iter, 4) = toc
    
end