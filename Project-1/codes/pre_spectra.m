clear all;
close all;

acq_freq = 60000; %sampling frequency
dt = 1./acq_freq; % time interval between successive data points
lags_t = 1; % The window over which we want to look at the correlation. This is in seconds.
nu = 1.5e-5;
Lint = 1800;
Nscales = 50;


fn = ['../flow1/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

[Pxx,F] = pwelch(un,Nscales.*Lint,1,acq_freq);
k = 2.*pi.*F/um;
diss_s = k.*k.*Pxx;

figure(1);
semilogx(k,k.*Pxx,'b-');
hold on;
diss = 15.*nu.*trapz(k,diss_s)

pause;
fn = ['../flow2/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

[Pxx,F] = pwelch(un,Nscales.*Lint,1,acq_freq);

figure(2);
k = 2.*pi.*F/um;
diss_s = k.*k.*Pxx;
semilogx(k,k.*Pxx,'r-');
diss = 15.*nu.*trapz(k,diss_s)


