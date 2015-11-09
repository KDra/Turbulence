%% example statistical analysis
% first look at the data
% t.g.thomas@soton.ac.uk (2015)

% filename containing the data
fn = ['../Project-1/flow2/u1_pos_11_burst1.bin']; 

% open the file, binary, and read it
fid = fopen(fn,'rb');               % rb=binary
u = fread(fid,inf,'float');         % read as floats
n = length(u);                      % number of samples

% take a look at the time series
dt = 1/60000;                       % sample interval [s]
T = (1./60000)*n;                   % sampling period [s]
t = (0:n-1)*dt;

figure(1)
hold off
plot(t,u);
xrange = [1, 1.1];
xlim(xrange);                       % zoom-in 
xlabel('t');
ylabel('u(t)');

% histogram and PDF
uvals = [0:0.05:20];                % grid of u-values
counts = hist(u,uvals);             % count into bins
area = trapz(uvals,counts);         % area under curve
pdf = counts./area;                 % normalise to recover PDF

% plots
figure(2)
hold off
plot(uvals,pdf)
xlabel('u')'
ylabel('PDF');  
hold on

% mean and variance (central moments)
U = mean(u)
U2 = var(u)
sigma = std(u)
S = skewness(u)
K = kurtosis(u)

% calc. central moments from the PDF
% should match those above
trapz(uvals,uvals.*pdf)
trapz(uvals,(uvals-U).^2.*pdf)
trapz(uvals,(uvals-U).^3.*pdf)/sigma^3
trapz(uvals,(uvals-U).^4.*pdf)/sigma^4

% normal PDF for comparison, same mean and std-dev
npdf = (1/sqrt(2*pi)/sigma)*exp(-(uvals-U).^2/2/sigma^2);
plot(uvals,npdf);
