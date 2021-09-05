%% MOC for reflected expansion wave in shock-tube
% This script plots the x-t graph for the reflection of the expansion wave
% in a 1-D shock-tube using method of characteristics (MOC).

%% Initialisation
% This section defines the properties at the driving the driven sections.
clear

n = 4;    % solver resolution
n_plot = 4;   % number of lines to plot

gamma = 1.4;    % adiabatic constant of air
R = 287;    % gas constant of air

PR41 = 5;    % pressure ratio = P4/P1
T1 = 300;   % driving section temperature
T4 = 300;    % driven section temperature
l = 3;    % length of driving section

%% Shock tube relations
% This section uses the shocktube relations to evaluate flow quantities in
% regions 3 and 4.

a1 = (gamma*R*T1)^0.5;
a4 = (gamma*R*T4)^0.5;

% iteratively solving for PR21
PR21 = PR41;
f = 1;
while abs(f)>0.1
    f = -PR41 + PR21*(1 - ((gamma-1)*(a1/a4)*(PR21-1)/...
        ((2*gamma*(2*gamma + (gamma+1)*(PR21-1)))^0.5)))^(-2*gamma/(gamma-1));
    if f>0
        PR21 = PR21 - 0.01;
    else
        PR21 = PR21 + 0.01;
    end
end

PR34 = PR21/PR41;

T2 = T1*PR21*((gamma+1)/(gamma-1) + PR21)/(1 + PR21*(gamma+1)/(gamma-1));
a2 = (gamma*R*T2)^0.5;
T3 = T4*(PR34^((gamma-1)/gamma));
a3 = (gamma*R*T3)^0.5;
u3 = 2*(a4 - a3)/(gamma-1);

%% MOC solver routine
% This section uses the method-of-characteristics to evaluate behaviour of
% the reflected expansion wave in the x-t plane.

xbyt = -a4:((u3 - a3 + a4)/(n-1)):(u3 - a3);    % equally spacing (x/t) between (-a4) and (u3-a3)
jn = (a4 + xbyt)*4/(gamma+1) - 2*a4/(gamma-1);  % vector containing j-minus invariants
jp = -jn;   % vector containing j-plus invariants
u = zeros(n);   % (n X n) matrix containing local u at points of intersection
a = u; x = u; t = u; mp = u; mn = u;    % local a, c-minus slope, c-plus slope and (x,t) coordinates

for i = 1:n
    for j = i:n
        u(i,j) = (jp(i) + jn(j))/2;
        a(i,j) = (jp(i) - jn(j))*(gamma-1)/4;
        mp(i,j) = 1/(u(i,j) + a(i,j));
        mn(i,j) = 1/(u(i,j) - a(i,j));
        
        if i==1 % for points lying on 1st reflected line
            x2 = 0; t2 = 0; m2 = mn(i,j);
        else
            x2 = x(i-1,j); t2 = t(i-1,j);
            m2 = (mn(i-1,j) + mn(i,j))/2;
        end
        
        if j==i % for points lying at wall
            x(i,j) = -l;
            t(i,j) = m2*(x(i,j) - x2) + t2; 
        else
            x1 = x(i,j-1); t1 = t(i,j-1);
            m1 = (mp(i,j-1) + mp(i,j))/2;
            
            x(i,j) = (m1*x1 - m2*x2 + t2 - t1)/(m1 - m2);
            t(i,j) = m1*(x(i,j) - x1) + t1;
        end
    end
end

%% Plots
% This section generates the plots of the characteristic lines on the x-t
% plane.
cx = zeros(1,n+2); ct = cx;     % characteristic lines to plot in x-t plane

for k = [(1:round((n-1)/(n_plot-1)):(n-1)),n]
    cx(1) = 0; ct(1) = 0;
    for j = 1:k
        cx(j+1) = x(j,k);
        ct(j+1) = t(j,k);
    end
    for j = (k+1):n
        cx(j+1) = x(k,j);
        ct(j+1) = t(k,j);
    end
    cx(n+2) = 0;
    ct(n+2) = -mp(k,n)*cx(n+1) + ct(n+1);
    plot(cx,ct,'LineWidth',2,'Color','k');
    xlabel('x'); ylabel('y'); title('Expansion Wave on xt Plane')
    hold on
end     

% end of script
