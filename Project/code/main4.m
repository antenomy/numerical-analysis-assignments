
%% Question 1

TV_x = linspace(0.6, 1);
TV_y = 0.6;

OMEGA = 30;

S0 = @(r) cos(24*sqrt(r)).*exp(-900*r);
S = @(amplitude, x, y, source_x, source_y) amplitude * S0((x-source_x).^2+(y-source_y).^2);

y_shift = 0.6; 

res_array = move_source(S, OMEGA, 200, TV_x, TV_y);
minimum = min(res_array, [], "all");
disp(minimum);
[min_x, min_y] = find(res_array==minimum);
S_min = @(x, y) S(1, x, y, TV_x(min_x), TV_y);
        
[bound, sol] = hhsolver(OMEGA, S_min, 200);

plotFields(bound, sol, S_min)

plotSoundRatio_x(res_array, TV_x)




function plotSoundRatio_x(A, x_shift) % A and x_shift are equally sized arrays
    
    plot(x_shift, A, 'b-')
    hold on;
    grid on;
    hold off;
end




%% Question 2

Q2_X_SHIFT = linspace(0.6, 1, 5); %= linspace(0.66, 0.7);
Q2_Y_SHIFT = 0.6;

OMEGA = 30;

S0 = @(r) cos(24*sqrt(r)).*exp(-900*r);
S = @(amplitude, x, y, source_x, source_y) amplitude * S0((x-source_x).^2+(y-source_y).^2);

res_array = move_source(S, OMEGA, 200, Q2_X_SHIFT, Q2_Y_SHIFT);

minimum = min(res_array, [], "all");
[min_x, min_y] = find(res_array==minimum);
a = Q2_X_SHIFT(min_x-1);
b = Q2_X_SHIFT(min_x+1); % We just assume these exist

disp(minimum)
disp(min_x)

disp(min_y)
disp(a)
disp(b)

TOL = 10e-4;

function result = func1(x, S, OMEGA, Q)
    res_array = move_source(S, OMEGA, 1000, x, y);
    result = res_array(1,1);
end

funcWrapper = @(x) func1(x, S);

final_x = goldenSectionSearch(@funcWrapper, a, b, TOL);

disp(minimum)
disp(final_x)
disp(funcWrapper(final_x))

% Golden section search
function result = goldenSectionSearch(func, a, b, tolerance)
    while (a-b) > tolerance
        g = (sqrt(5) - 1) / 2;
        if func(a + (1 - g)*(b-a)) < func(a + g*(b-a))
            b = a + g*(b-a);
        else
            a = a + (1-g) * (b-a);
        end
    end
    result = (a+b)/2;
end



%% Question 3




%%%%%%% Helper functions

function return_array = move_source(S, omega, N, x_shift, y_shift) % S Function, x_shift and y_shift arrays
    return_array = zeros(length(x_shift), length(y_shift));
    for i_x=1:length(x_shift)
        for i_y=1:length(y_shift)
            S_shifted = @(x, y) S(1, x, y, x_shift(i_x), y_shift(i_y));
            
            [bound, sol] = hhsolver(omega, S_shifted, N);
            w = find(sol.x<=0.25 & sol.y>=0.5);
            A = max(abs(sol.u(w)))/max(abs(sol.u(:)));

            return_array(i_x,i_y) = A;
            
            
        end
    end
end

function plotFields(BB, Sol, sourceModel)
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
    contour(Sol.x,Sol.y,Sol.u,20)
    plot3(BB.x,BB.y,zeros(size(BB.x)),'k-','LineWidth',2)
    hold off
        
    figure;
    mesh(Sol.x,Sol.y,Sol.u)
end
