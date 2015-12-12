clear all;
close all;

acq_freq = 60000; %sampling frequency
dt = 1./acq_freq; % time interval between successive data points
lags_t = 1; % The window over which we want to look at the correlation. This is in seconds.
nu = 1.5e-5;
Lint = 400;
Nscales = 50;
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
um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

[Pxx,F] = pwelch(un,Nscales.*Lint,1,acq_freq);

figure(1);
loglog(F,Pxx,'b-');
hold on;
% pause;

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
um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

[Pxx,F] = pwelch(un,Nscales.*Lint,1,acq_freq);

figure(1);
loglog(F,Pxx,'r-');


