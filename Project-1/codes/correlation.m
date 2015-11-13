clear all;
close all;

acq_freq = 60000; %sampling frequency
dt = 1./acq_freq; % time interval between successive data points
lags_t = 0.2; % The window over which we want to look at the correlation. This is in seconds.
Nfiles = 5;
ucor = [];
% loop over ensemble files
for i = 1:Nfiles

    % open the file, binary, and read samples 
    fn = sprintf('../flow1/u1_pos_11_burst%d.bin', i);
    %fn = sprintf('../flow2/u1_pos_11_burst%d.bin', i);
    fid = fopen(fn,'rb'); 
    u = fread(fid,inf,'float'); 
    %fprintf(1,'Read %d samples from file %s\n', n, fn);
    ucor = [ucor; u];   
end

u = ucor;
len = length(u); %finding the number of points in the vector
total_time = len.*dt; %establishing the total length (or time) of the signal (since we know acquisition frequency)

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal
us = std(un);

lags_n = floor(lags_t.*acq_freq); %calculating the number of points over which the correlation should be calculated

%[c,lags] = xcorr(un,un,'unbiased'); %calcuates the autocorrelation
[c,lags] = autocorr(un, lags_n);
rho = c;%./(us.*us); %normalises the autocorrelation

plot(lags.*dt,rho,'b-');
hold on;
%
f = find(lags>=0);
dr = [fliplr(rho(2:3)); rho(1:3)];
der = [-1/12 	4/3 	-5/2 	4/3 	-1/12];
tl1 = sqrt(-2 * dt^2 / sum(dr .* der'))
Lf1 = trapz(lags,rho)%/acq_freq

%%
ucor = [];
% loop over ensemble files
for i = 1:Nfiles

    % open the file, binary, and read samples 
    %fn = sprintf('./flow1/u1_pos_11_burst%d.bin', i);
    fn = sprintf('../flow2/u1_pos_11_burst%d.bin', i);
    fid = fopen(fn,'rb'); 
    u = fread(fid,inf,'float'); 
    %fprintf(1,'Read %d samples from file %s\n', n, fn);
    ucor = [ucor; u];   
end

u = ucor;
len = length(u); %finding the number of points in the vector
total_time = len.*dt; %establishing the total length (or time) of the signal (since we know acquisition frequency)

um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal
us = std(un);

lags_n = floor(lags_t.*acq_freq); %calculating the number of points over which the correlation should be calculated

%[c,lags] = xcorr(un,lags_n,'unbiased'); %calcuates the autocorrelation
[c,lags] = autocorr(un, lags_n);
rho = c;%./(us.*us); %normalises the autocorrelation


plot(lags.*dt,rho,'r-');

dr = [fliplr(rho(2:3)); rho(1:3)];
tl2 = sqrt(-2 * dt^2 / sum(dr .* der'))
Lf2 = trapz(lags,rho)% /acq_freq