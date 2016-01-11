%% 1D Convection Diffusion Equation
% t.g.thomas@soton 23/Nov/2015

function [flag] = wave(c, nu)
flag = 1;
% debug
% fprintf(1,'Start: wave.m\n');

% numerical parametes
N = 2048;                               % number of grid points
X0 = -5.0;                              % start of domain
X1 = 10.0;                              % end of domain
L = X1 - X0;                            % length of domain
%
dt = 0.003                            % time interval
dx = L/N                               % mesh interval

% physical parameters
%c = 1.0;                                % convection velocity
%nu = 1.0e-5;                            % diffusion

% artificial viscosity (Lax-Wendroff)
%nu = 0.501*c^2*dt;                       % Euler-CDif2/CDif4
%nu = 0.03*c^4*dt^3/dx^2;                 % RK2-CDif2/CDif4

% solution
t = 0;                                  % time
f = zeros(1,N);                         % solution

% grid
x = X0 + (0:N-1)*dx;                    % uniform grid

% initial solution (at t=0)
sigma = 0.05;                           % gaussian semi-width
f = gaussian(0, sigma, x);

% integral constants
% fprintf(1,'Initial, Area = %5.3e\n', dx*sum(f));
% fprintf(1,'Initial, Energy = %5.3e\n', dx*sum(f.^2)/2);
miE = dx*sum(f.^2)/2;
% Courant numbers
% fprintf(1,'Convection CFL: c*dt/dx = %5.3e\n', c*dt/dx);
% fprintf(1,'Diffusion CFL: nu*dt/dx^2 = %5.3e\n', nu*dt/dx^2);
f1=f;
%  time steps
Nstep = 2500;                             % number of time steps

% record time-series data
T = zeros(1,Nstep);
A = zeros(1,Nstep);
E = zeros(1,Nstep);

for k = 1:Nstep
    % Edit: select scheme below
	% Change time-steps method
    %f = tstep_Euler(@dfdt_diff2, f, c, nu, dt, dx, N); % unstable
    %f = tstep_Euler(@dfdt_diff4, f, c, nu, dt, dx, N); 
    f = tstep_RK2(@dfdt_diff2, f, c, nu, dt, dx, N);
    %f = tstep_RK2(@dfdt_diff4, f, c, nu, dt, dx, N);
    
    % logging
    t = t + dt;
    A(k) = dx*sum(f);
    E(k) = dx*sum(f.^2)/2;
    T(k) = t;
end

% evaluate the solution
% integral 'constants'
% fprintf(1,'Area = %5.3e\n', dx*sum(f));
% fprintf(1,'Energy = %5.3e\n', dx*sum(f.^2)/2);

en = dx*sum(f.^2)/2 - miE;
% dedt = (E(2:end)-E(1:end-1))/dt;
% dedt2 = (dedt(2:end) - dedt(1:end-1))/dt;
% df = dfdt_diff4(f, c, nu, dx, N);
% df2 = dfdt_diff4(df, c, nu, dx, N); 
if en>0 || isinf(en) || isnan(en) %sum(abs((f1-f)))>10
    flag = 0;
end
% % plots
% figure(1)
% hold off
% plot(x,f)
% hold on
% plot(x,gaussian(0+c*T(Nstep),sigma,x),'--r');
% xlabel('x');
% ylabel('f');
% 
% % repeat + zoom in
% figure(2)
% hold off
% plot(x,f)
% hold on
% plot(x,gaussian(0+c*T(Nstep),sigma,x),'--r');
% xlim([c*T(Nstep)-1,c*T(Nstep)+1]);
% xlabel('x');
% ylabel('f');
% 
% figure(3)
% hold off
% plot(T,A);
% hold on
% plot(T,E,'r');
% xlabel('t');
% ylabel('A(t), E(t)');
% hold off
end                                     % end of wave()



% tstep_Euler - advance the solution over a time step using Euler scheme
function [f1] = tstep_Euler(dfdt, f0, c, nu, dt, dx, N)

f1 = zeros(1,N);
f1 = f0 + dt*dfdt(f0, c, nu, dx, N);
end % tstep_Euler 

% tstep_RK2 - advance the solution over a time step using RK2
function [f1] = tstep_RK2(dfdt, f0, c, nu, dt, dx, N)

ftmp = zeros(1,N);
f1 = zeros(1,N);

% predictor
ftmp = f0 + 0.5*dt*dfdt(f0, c, nu, dx, N);
% corrector
f1 = f0 + dt*dfdt(ftmp, c, nu, dx, N);
end % tstep_RK2


% dfdt_diff2 - return df/dt calculated from the convection-diffusion equation
function [df] = dfdt_diff2(f, c, nu, dx, N)

df = zeros(size(f));

% interior points
for j = 2:N-1
    df(j) =  -(0.5*c/dx)*(f(j+1) - f(j-1)) ...
        + (nu/dx^2)*(f(j+1) - 2*f(j) + f(j-1));
end

% wrap around boundary
df(1) = -(0.5*c/dx)*(f(2) - f(N)) ...
    + (nu/dx^2)*(f(2) - 2*f(1) + f(N));
df(N) = -(0.5*c/dx)*(f(1) - f(N-1)) ...
    + (nu/dx^2)*(f(1) - 2*f(N) + f(N-1));
end % dfdt_diff2


% dfdt_diff4 - return df/dt calculated from the convection-diffusion equation
function [df] = dfdt_diff4(f, c, nu, dx, N)

df = zeros(size(f));

% interior points
for j = 3:N-2
    df(j) =  -(c/(12*dx))*(-f(j+2) + 8*f(j+1)- 8*f(j-1) + f(j-2)) ...
        + (nu/(12*dx^2))*(-f(j+2) + 16*f(j+1) - 30*f(j) + 16*f(j-1) - f(j-2));
end

% wrap around boundary
df(1) = -(c/(12*dx))*(-f(3) + 8*f(2) - 8*f(N) + f(N-1)) ...
    + (nu/(12*dx^2))*(-f(3) + 16*f(2) - 30*f(1) + 16*f(N) - f(N-1));
df(2) = -(c/(12*dx))*( -f(4) + 8*f(3) - 8*f(1) + f(N)) ...
    + (nu/(12*dx^2))*(-f(4) + 16*f(3) - 30*f(2) + 16*f(1) - f(N));
df(N) = -(c/(12*dx))*(-f(2) + 8*f(1) - 8*f(N-1) + f(N-2)) ...
    + (nu/(12*dx^2))*(-f(2) + 16*f(1) - 30*f(N) + 16*f(N-1) - f(N-2));
df(N-1) = -(c/(12*dx))*(-f(1) + 8*f(N) - 8*f(N-2) + f(N-3)) ...
    + (nu/(12*dx^2))*(-f(1) + 16*f(N) - 30*f(N-1) + 16*f(N-2) - f(N-3));
end % dfdt_diff4


% gaussian - return a Gaussian profile of semi-width sigma located at x=xc
function [f] = gaussian(xc, sigma, x)
% xc = location x of the peak 
% sigma = Gaussian semi-width
% x = {x1, x2, ..., xn} locations to return the profile
% f = {f1, f2, ..., fn} values of f at x
f = zeros(size(x));

for j = 1:length(x)
  f(j) = exp(-(x(j)-xc).^2/(2*sigma.^2))/sqrt(2*pi)/sigma;
end
end % gaussian


% random - return