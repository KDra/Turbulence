%% example time-series analysis
% t.g.thomas@soton.ac.uk (Oct/2015)

% filename containing the data
fn = ['../flow2/u1_pos_11_burst1.bin']; 

% open the file, binary, and read it
fid = fopen(fn,'rb');               % rb=binary
u = fread(fid,inf,'float');         % read as floats
n = length(u);                      % number of samples

% take a look at the time series
SR = 60000;                         % sample rate [S/s]
dt = 1/SR;                          % sample interval [s]
T = dt*n;                           % sampling period [s]
t = (0:n-1)*dt;

% figure(1)
% hold off
% plot(t,u);
% xrange = [1, 1.1];
% xlim(xrange);                       % zoom-in 
% xlabel('t');
% ylabel('u(t)');

% statistics
uave = mean(u);                     % average velocity
sigma = std(u);                     
u = u - uave;                       % remove mean

% Fourier PSP - square window
df = 1/T;                           % spectral bin size
f = (0:n-1)*df;                     % spectral bins
F = fft(u)/n;                       % FFT, normaise power                    
S = F.*conj(F)/(df);                % power density per Hz
n2 = floor(n/2);                    % don't plot beyond n/2

% sanity check - power in t-domain and f-domain
var(u)                              % signal power in time domain
trapz(f,S)                          % signal power in spectral domain

figure(2)
loglog(f(1:n2),2*S(1:n2));        % one-sided spectrum, so 2x
hold on
plot(f(1:n2),smooth(2*S(1:n2),7));  % smooth, average adjacent bins

% Fourier PSP - Hamming window, one-side,Welsh averaging, overlapping windows
windowSize = 6.0*SR;                % window length, [n samples]
overlapSize = 3.0*SR;               % overlap length [n samples]
[Pxx,F] = pwelch(u, windowSize, overlapSize, 4096, SR);

figure(3);
%loglog(F,F.*Pxx,'r-');             % pre-mult spectrum
loglog(F,Pxx,'r-');


