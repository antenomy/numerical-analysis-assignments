%% Interpolation
h = 0.25; % Får ej ändras i koden nedan
%% Linjär interpolation
%kraven uppfyllda
[t,x,y,vx,vy] = kastbana(h);
[xmax_linjar, ymax_linjar, X_nedslag_linjar] = styckvis_linjar_interpolation(x,y)



plot(xmax_linjar, ymax_linjar, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
 plot(X_nedslag_linjar, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
%% Kvadratisk interpolation
%kraven uppfyllda
[t,x,y,vx,vy] = kastbana(h);
[xmax_kvad, ymax_kvad, X_nedslag_kvad] = styckvis_kvadratisk_interpolation(x,y)

 plot(xmax_kvad, ymax_kvad, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(X_nedslag_kvad, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'r');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x,y,vx,vy] = kastbana(h)
%KASTBANA(H) beräknar banan för ett kast med en liten boll.
%
% Dynamiken ges av en ODE som inkluderar effekten av luftmotståndet,
%
% r'' = -g*ez-sigma*r'*|r'|/m.
%
% Funktionen beräknar bollens position och hastighet vid
% tidpunkter separerade med en given steglängd. Bollen kastas från
% (X,Y)=(0,2) med hastigheten 30 m/s i 45 graders vinkel uppåt.
%
% Syntax:
%
% [T,X,Y,VX,VY] = KASTBANA(H)
%
% H - Steglängden mellan tidpunkterna.
% T - Vektor med tidpunkter där bollens position och hastighet beräknats.
% X, Y - Vektorer med bollens x- och y-koordinater vid tidpunkterna.
% VX, VY - Vektorer med bollens hastigheter i x- och y-led vid tidpunkterna.



%% Tennisboll, specifikationer
m = 56e-3; % Massan (kg) = 56 gram
ra = 6.6e-2/2; % 6.6 cm in diameter
g=9.81; % Tyngdaccelerationen (m/s^2)
rho=1.2; % Luftens densitet (kg/m^3)
A=ra^2*pi; % Kroppens tvärsnittsarea (m^2)
Cd=0.47; % Luftmotståndskoefficient,
% "drag coefficient" (dimensionslös)
% Läs mer på http://en.wikipedia.org/wiki/Drag_coefficient
sigma = rho*A*Cd/2; % Totala luftmotståndet
T = 5; % Sluttid
v0 = 32; % Utkasthastighet
al = pi/4; % Utkastvinkel
% Begynnelsevärden
r0 = [0 2]'; % Position
r1 = [v0*cos(al) v0*sin(al)]'; % Hastighet
% ODEns högerled
f = @(u) [u(3:4); -u(3:4)*norm(u(3:4),2)*sigma/m - [0;g]]; % RHS
u = [r0;r1];
U = u';
t = 0:h:T;
% Runge-Kutta 4
for tn=t(1:end-1)
s1 = f(u);
s2 = f(u + h/2*s1);
s3 = f(u + h/2*s2);
s4 = f(u + h*s3);
u = u + h/6*(s1 + 2*s2 + 2*s3 + s4);
U = [U; u'];
end
x = U(:,1);
y = U(:,2);
vx = U(:,3);
vy = U(:,4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_max, y_max, x_nedslag] = styckvis_linjar_interpolation(x,y)
figure(1)
    for i=1:length(x)-1
       x0 = x(i);x1 = x(i+1);y0 = y(i);y1 = y(i+1);
       plot([x0,x1], [y0,y1], '-r', 'LineWidth', 2)
        hold on; 

        if y(i) * y(i+1) <= 0
            k = (y1 - y0)/(x1 - x0);
            
            f = @(X) X - ((k.*X + y0 - k.*x0)/(k)); %newtons konvergens formel
               
            x_nedslag = 48.911; %startgissning från plot
            tau = 10^-10;
            fel = abs(f(x_nedslag) - x_nedslag);
            format long;
            while fel >= tau
                X = f(x_nedslag);%sätter nya värdet för nästa itteration
                fel = abs(f(X)-X); %uppdaterar felet
                
                x_nedslag = X; %uppdaterar, för att få rätt värden till nästa itteration
            end
        end
    end
    [y_max, y_max_index] = max(y);
     
    x_max = x(y_max_index);
   
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_max, y_max, x_nedslag] = styckvis_kvadratisk_interpolation(x,y)
figure(2)
    n = length(x);
       
    V = ones(n,3);
    V(: , 2) = x(:);
    V(:, 3) = x(:).^2;
    
    C = V\y;
    
    c0 = C(1);
    c1 = C(2);
    c2 = C(3);
    x = linspace(0, 70, 100);
    f = @(x) c0 + c1.*x + c2.*x.^2;
    y = f(x);
    plot(x, y, 'b-', 'LineWidth', 2);
    hold on
    

    f = @(X) X - ((c0 + c1.*X + c2.*X.^2)/(c1 + 2.*X.*c2));
    x_nedslag = 48.911; %startgissning från plot
    tau = 10^-10;
    fel = abs(f(x_nedslag) - x_nedslag);
    format long;
    while fel >= tau
        X = f(x_nedslag);%sätter nya värdet för nästa itteration
        fel = abs(f(X)-X); %uppdaterar felet
                    
        x_nedslag = X; %uppdaterar, för att få rätt värden till nästa itteration       
    end



    x_max = 24.0404; %startgissning från plot
    f = @(X) X - ((c1 + 2.*X.*c2)/(2.*c2)); %Newton Rapshon metoden
    while fel >= tau 
        X = f(x_max);
        fel = abs(f(X)-X);

        x_max = X;
    
    end

    f = @(x) c0 + c1.*x + c2.*x.^2;
    y_max = f(x_max);
    
   

end