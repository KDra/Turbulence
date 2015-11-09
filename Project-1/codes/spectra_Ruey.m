clear all;
close all;

acq_freq = 10000; %sampling frequency (I think in your case, this is 10 kHz)
dt = 1./acq_freq; % time interval between successive data points
win_size = 20000; %This is the number of points over which the spectra is calculated. I have used 20000 points. This means we break up the signal in to lots of 2 second blocks and see what the spectra shows. 

Union = 10; %This is the freestream velocity and is equal to 10 m/s
D = 0.05; % This is the diameter or size of the cylinder (0.05m). 


fn = ['../flow1/u1_pos_11_burst1.bin']; %assigning a file name to read
fid = fopen(fn,'rb'); % opening a file in binary form so that you can read the file
u = fread(fid,inf,'float'); %reading the data in the file in to a vector - this reads all the data in to the this vector


um = mean(u); %calculate the mean of the signal
un = u-um;%calculate the fluctuation of the signal

[Pxx,F] = pwelch(un,win_size,1,acq_freq);


figure(1);
semilogx(F.*D/Uinf,F.*Pxx./Uinf.^2,’b-‘);


