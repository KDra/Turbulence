%% example time-series analysis
% t.g.thomas@soton.ac.uk (Oct/2015)

% filename containing the data
fn = ['../flow2/u1_pos_11_burst1.bin']; 

% open the file, binary, and read it
fid = fopen(fn,'rb');               % rb=binary
u = fread(fid,inf,'float');         % read as floats
n = length(u);                      % number of samples

% take a look at the time series
SR = 60000.0;                       % sample rate [S/s]
dt = 1/SR;                       % sample interval [s]
T = (1./SR)*n;                   % sampling period [s]
t = (0:n-1)*dt;

figure(1)
hold off
plot(t,u);
xrange = [1, 1.1];
xlim(xrange);                       % zoom-in 
xlabel('t');
ylabel('u(t)');

% statistics
uave = mean(u);                     % average velocity
sigma = std(u);                     
u = u - uave;                       % remove mean

% correlations
maxTlag = 1.0                       % restrict maximum lag [s]
maxnlag = floor(maxTlag*SR); 
[c,lags] = xcorr(u,maxnlag,'unbiased'); % unbiased=correct reduced overlap
R = c./sigma^2;                     % normalise

% integral scales, Taylor frozen turbulence idea to convert to length 
integralT = trapz(lags*dt, R)/2          % one-sided area
fprintf(1,'Integral Time Scale = %f s\n', integralT);
fprintf(1,'Integral Length Scale = %f m\n', integralT*uave);


% Taylor microscale
izero = maxnlag + 1;                % index of zero lag, R(izero)=1
d2Rdt2 = (R(izero+1) -2*R(izero) + R(izero-1))/(dt*dt);
taylorT = sqrt(-2.0/d2Rdt2);
fprintf(1,'Taylor microscale (time) = %f s\n', taylorT);
fprintf(1,'Taylor microscale (length) = %f m\n', taylorT*uave);

figure(2)
hold off
plot(lags*dt,R);
hold on
%plot(lags*dt, 1-(dt^2*lags.^2/taylorT^2),'r');
ylim([-0.1,1.1]);
xlim([-0.1,0.1]);

figure(3)
hold off
plot(lags*dt,R);
hold on
plot(lags*dt, 1-(dt^2*lags.^2/taylorT^2),'r');
ylim([-0.1,1.1]);
xlim([-0.001,0.001]);


