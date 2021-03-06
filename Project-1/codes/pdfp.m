clear all;
close all;

acq_freq = 60000; %sampling frequency
dt = 1./acq_freq; % time interval between successive data points
lags_t = 1; % The window over which we want to look at the correlation. This is in seconds.



fn = ['../flow1/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

ax = [-6:.01:6];
ay = hist(un,ax);

ay = ay./trapz(ax,ay);

plot(ax,ay,'b-');
hold on; 

f1_Moment1 = um
f1_Moment2 = var(un)
f1_Moment3 = skewness(un)
f1_Moment4 = kurtosis(un)

pause;

fn = ['../flow2/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

ax = [-6:.01:6];
ay = hist(un,ax);

ay = ay./trapz(ax,ay);

plot(ax,ay,'r-');

f2_Moment1 = um
f2_Moment2 = var(un)
f2_Moment3 = skewness(un)
f2_Moment4 = kurtosis(un)

