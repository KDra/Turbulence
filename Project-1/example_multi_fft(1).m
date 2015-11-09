%% example time-series analysis with ensemble average
% t.g.thomas@soton.ac.uk (Nov/2015)

clear all

% Example of how to ensemble average power spectra
% over multiple data files.

nfft = 4096;                        % length of FFT used in pwelch
SR = 60000;                         % sample rate [S/s]
dt = 1/SR;                          % sample interval [s]
Twindow = 1.0;                      % FFT window interval [S]
Nfiles = 5;                         % number of files

% accumulate PSD from each file into these Fourier bins
big_Pxx = zeros(1+nfft/2,1);        % pwelch returns array of this size

% loop over ensemble files
for i = 1:Nfiles
   
    % read samples from file
    %fn = sprintf('../flow1/u1_pos_11_burst%d.bin', i);
    fn = sprintf('../flow2/u1_pos_11_burst%d.bin', i);
    fid = fopen(fn,'rb');           % rb=binary
    u = fread(fid,inf,'float');     % read as floats
    n = length(u);         
    fprintf(1,'Read %d samples from file %s\n', n, fn);
    
    % remove mean
    u = u - mean(u);
    % check variance in t-domain
    fprintf(1, 'Variance %f\n', var(u));
    
    % PSD, sliding window with 50% overap
    windowSize = Twindow*SR;        % window length, [n samples]
    overlapSize = 0.5*windowSize;   % overlap length [n samples]
    [Pxx,F] = pwelch(u, windowSize, overlapSize, nfft, SR);
    
    % accumulate PSD
    big_Pxx = big_Pxx + Pxx;
end
big_Pxx = big_Pxx/Nfiles;

% check power in Fourier domain
fprintf(1, 'Combined Fourier Power %f\n', trapz(F,big_Pxx));

figure(10);
%loglog(F,F.*big_Pxx,'r-');         % pre-mult spectrum
loglog(F,(F).*big_Pxx,'b-');             % normal energy spectrum


