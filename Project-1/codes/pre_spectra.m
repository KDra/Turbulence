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
k = 2.*pi.*F/um;
diss_s = k.*k.*Pxx;
diss1 = 15.*nu.*trapz(k,k.*Pxx)
lambda = sqrt(2*trapz(F, Pxx)/trapz(k,k.*Pxx))
%(nu^3/diss1)^0.25
L1 = k(find(max(Pxx.*k)==Pxx.*k))
%figure(1);
loglog(k,k.*Pxx,'b-');
hold on;

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

%figure(2);
k = 2.*pi.*F/um;
diss_s = k.*k.*Pxx;
loglog(k,k.*Pxx,'r-');
diss2 = 15.*nu.*trapz(k,diss_s/dt)
lambda = sqrt(2*trapz(F, Pxx)/trapz(k,k.*Pxx))
%(nu^3/diss2)^0.25

L2 = k(find(max(Pxx.*k)==Pxx.*k)) %Pxx(1)/4/um^2
%