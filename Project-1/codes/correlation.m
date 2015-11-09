clear all;
close all;

acq_freq = 60000; %sampling frequency
dt = 1./acq_freq; % time interval between successive data points
lags_t = 1; % The window over which we want to look at the correlation. This is in seconds.



fn = ['../flow1/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector
len = length(u); %finding the number of points in the vector
total_time = len.*dt; %establishing the total length (or time) of the signal (since we know acquisition frequency)

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal
us = std(un);

lags_n = floor(lags_t.*acq_freq); %calculating the number of points over which the correlation should be calculated

[c,lags] = xcorr(un,lags_n,'unbiased'); %calcuates the autocorrelation

figure(1);
subplot(1,2,1);
plot(lags.*dt,c,'b-'); %plots the autocorrelation over a given window
hold on;

rho = c./(us.*us); %normalises the autocorrelation

subplot(1,2,2);
plot(lags.*dt,rho,'b-');
hold on;
% 
Lf1 = trapz(lags,rho)
pause

fn = ['../flow2/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector
len = length(u); %finding the number of points in the vector
total_time = len.*dt; %establishing the total length (or time) of the signal (since we know acquisition frequency)

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal
us = std(un);

lags_n = floor(lags_t.*acq_freq); %calculating the number of points over which the correlation should be calculated

[c,lags] = xcorr(un,lags_n,'unbiased'); %calcuates the autocorrelation

figure(1);
subplot(1,2,1);
plot(lags.*dt,c,'r-'); %plots the autocorrelation over a given window


rho = c./(us.*us); %normalises the autocorrelation
Lf2 = trapz(lags,rho)


figure(1);
subplot(1,2,2);
plot(lags.*dt,rho,'r-');
% 


