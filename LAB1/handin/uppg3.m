%% Assignment a)

load('eiffel4.mat');
load('eiffel1.mat');

trussplot(xnod,ynod,bars, 'w')
hold on

matrix_size = height(A);
b = zeros(matrix_size,1); 
bel_index = round(matrix_size/3);
b(bel_index) = 1;

x = A \ b;

xbel = xnod + x(1:2:end); 
ybel = ynod + x(2:2:end);

trussplot(xbel,ybel,bars, 'r'); 
hold on
plot(xnod(bel_index), ynod(bel_index), 'y*', 'MarkerSize', 15);

hold off


%% Assignment b)

iter = 1;
unkown_variables = zeros(1, 4);
compute_times = zeros(1, 4);

for file = {'eiffel1.mat', 'eiffel2.mat', 'eiffel3.mat', 'eiffel4.mat'}
    disp(file{1})
    load(file{1})
    matrix_size = height(A);
    unkown_variables(iter) = matrix_size;

    b = randn(matrix_size,1);

    tic;
    x = A \ b;
    compute_times(iter) = toc;

    iter = iter + 1;
end

loglog(unkown_variables, compute_times, '-o')

xlabel('Unkown Variables');
ylabel('Compute Times (s)');
grid on

%% Assignment c)





%% Assignment d) - Time Table


% Skapa en 4x4-matris T som innehåller beräkningstiderna.
% Raderna ska motsvara de olika modellerna (eiffel1-eiffel4) och
% kolumnerna de olika metoderna, ordnade som "Naiv", "LU",
% "Gles" och "Gles LU".

% Följande kod skapar en snygg tabell med resultaten:

T = zeros(4, 4);

tab=array2table(T,'VariableNames',{'Naiv' 'LU' 'Gles' 'Gles LU'},'RowNames',{'eiffel1' 'eiffel2' 'eiffel3' 'eiffel4'});
disp(tab);

%% Two

x = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [jmin,jmax]=kanslighet(A,metod)


%  Indata:
%
%  A     - matrisen
%  metod - villken metod som används:
%          1 = Naiv metod
%          2 = LU-faktorisering
%
%  Utdata:
%
%  jmin - index för minst känsliga nod
%  jmax - index för mest känsliga nod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...

%end
