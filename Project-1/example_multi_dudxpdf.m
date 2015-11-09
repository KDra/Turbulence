%% example PDF with ensemble average
% t.g.thomas@soton.ac.uk (Nov/2015)

% Example of how to ensemble average PDFs
% over multiple data files.

clear all

uvals = [-1.0e4:50:1.0e4];                % grid of u-values
Nfiles = 5;                         % number of files
SR = 60000;                         % sample rate [S/s]
dt = 1/SR;                          % time interval [s]
Umean = 10.0362;                    % mean velocity
nu = 1.5e-5;                        % kinematic viscosity [m^2/s]

% accumulate counts from each file into these bins
big_counts = zeros(1,length(uvals));

% loop over ensemble files
for i = 1:Nfiles

    % open the file, binary, and read samples 
    %fn = sprintf('../flow1/u1_pos_11_burst%d.bin', i);
    fn = sprintf('../flow2/u1_pos_11_burst%d.bin', i);
    fid = fopen(fn,'rb'); 
    u = fread(fid,inf,'float'); 
    n = length(u);
    fprintf(1,'Read %d samples from file %s\n', n, fn);

    % derive dudx from u, use Taylor um
    for j = 2:n-1
        dudx(j) = 0.5*(u(j+1) - u(j-1))/(Umean*dt);
    end
    dudx(1) = (u(2) - u(1))/(Umean*dt);
    dudx(n) = (u(n) - u(n-1))/(Umean*dt);
    
    % print min/max dudx
    fprintf('dudx min%f nax %f\n', min(dudx), max(dudx)); 
    
    % histogram and PDF
    counts = hist(dudx,uvals);             % count into bins

    % accumulate counts
    big_counts = big_counts + counts;    
end
big_counts = big_counts/Nfiles;

% normalize for unit area
area = trapz(uvals,big_counts);         % area under curve
pdf = big_counts./area;                 % normalise to recover PDF

% plots
figure(2)
hold off
plot(uvals,pdf)
hold on

% calculate central moments from the PDF
U = trapz(uvals,uvals.*pdf)
variance = trapz(uvals,(uvals-U).^2.*pdf)
sigma = sqrt(variance)
skew = trapz(uvals,(uvals-U).^3.*pdf)/sigma^3
kurt = trapz(uvals,(uvals-U).^4.*pdf)/sigma^4

% dissipation
dissipation = 15.*nu.*variance

% plot fitted normal PDF
%npdf = (1/sqrt(2*pi)/sigma)*exp(-(uvals-U).^2/2/sigma^2);
%plot(uvals,npdf);