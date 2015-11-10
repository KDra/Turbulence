%% example PDF with ensemble average
% t.g.thomas@soton.ac.uk (Nov/2015)

%clear all

% Example of how to ensemble average PDFs
% over multiple data files.

% adjust bin limits to fit the data
uvals = [0:0.05:20];                % grid of u-values
Nfiles = 5;                         % number of files

% accumulate counts from each file into these bins
big_counts = zeros(1,length(uvals));

% loop over ensemble files
for i = 1:Nfiles

    % open the file, binary, and read samples 
    fn = sprintf('./flow1/u1_pos_11_burst%d.bin', i);
    %fn = sprintf('./flow2/u1_pos_11_burst%d.bin', i);
    fid = fopen(fn,'rb'); 
    u = fread(fid,inf,'float'); 
    n = length(u);
    fprintf(1,'Read %d samples from file %s\n', n, fn);

    % histogram and PDF
    counts = hist(u,uvals);             % count into bins

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

% plot fitted normal PDF
npdf = (1/sqrt(2*pi)/sigma)*exp(-(uvals-U).^2/2/sigma^2);
plot(uvals,npdf);
legend('Flow','Gaussian')