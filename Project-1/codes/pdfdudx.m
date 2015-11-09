clear all;
close all;

acq_freq = 60000; %sampling frequency
dt = 1./acq_freq; % time interval between successive data points
lags_t = 1; % The window over which we want to look at the correlation. This is in seconds.
nu = 1.5e-5;


fn = ['../flow1/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal
dudx = (diff(un(2:length(un)))-diff(un(1:length(un)-1)))./(2.*um.*dt);


ax = [min(dudx):10:max(dudx)];
ay = hist(dudx,ax);

ay = ay./trapz(ax,ay);

semilogy(ax,ay,'r-');
axis([-1e4 1e4 1e-7 1e-2]); 
hold on; 
f1_Moment1 = mean(dudx)
f1_Moment2 = var(dudx)
f1_Moment3 = skewness(dudx)
f1_Moment4 = kurtosis(dudx)

dissipation = 15.*nu.*f1_Moment2
eta = (nu.^3./dissipation).^(0.25)

pause;

fn = ['../flow2/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal
dudx = (diff(un(2:length(un)))+diff(un(1:length(un)-1)))./(2.*um.*dt);


ax = [min(dudx):10:max(dudx)];
ay = hist(dudx,ax);

ay = ay./trapz(ax,ay);

semilogy(ax,ay,'b-');

f2_Moment1 = mean(dudx);
f2_Moment2 = var(dudx)
f2_Moment3 = skewness(dudx)
f2_Moment4 = kurtosis(dudx)

dissipation = 15.*nu.*f2_Moment2

eta = (nu.^3./dissipation).^(0.25)



