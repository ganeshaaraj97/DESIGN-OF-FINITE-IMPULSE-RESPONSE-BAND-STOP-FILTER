%Generation of the required specification by the index number

I_A = 1; I_B = 8; I_C = 0;
%rad/s(below every parameter)
A_p = 0.03+(0.01*I_A); %max passband ripple
A_a = 45+I_B; %min stopband attenuation
lp = (I_C*100)+400; %lower passband edge
up = (I_C*100)+950; %upper passband edge
ls = (I_C*100)+500; %lower stopband edge
us = (I_C*100)+800; %upper stopband edge
sf = 2*((I_C*100)+1300); %sampling freqency
% A_p
% A_a
% O_p1
% O_p2
% O_a1
% O_a2
% O_s

%Obtaining the specifications of the filter

ltw = ls-lp; %lower transition width
utw = up-us; %upper transisiton width
ctw = min(ltw,utw); %critical transition width
low_cf = lp+ctw/2; %lower cutoff frequency
up_cf = up-ctw/2; %upper cutoff frequency
st = 2*pi/sf; %sampling period
% B_t1
% B_t2
% B_t
% O_c1
% O_c2
% T

%Obtaining the Kaiser Window

dp = ((10^(0.05*A_p)) -1)/(1+(10^(0.05*A_p)));
da = 10^(-0.05*A_a);
d = min(dp,da);
I_A = -20*log10(d);
Ap = 20*log10((1+d)/(1-d));
Aa = -20*log10(d);
if I_A<=21
 alpha = 0;
elseif I_A>21 && I_A<=50
 alpha = 0.5842.*(I_A-21).^0.4+0.07886.*(I_A-21);
else
 alpha = 0.1102.*(I_A-8.7);
end
%Calculating D
if I_A<=21
 D = 0.9222;
else
 D = (I_A-7.95)/14.36;
end
%Finding the order of the filter
N = ceil((sf*D/ctw)+1);
%Order of the filter should be odd
if mod(N,2) == 0
 N = N+1;
else
 N =N;
end
n = -(N-1)/2:1:(N-1)/2;
beta = alpha*sqrt(1-(2*n/(N-1)).^2);
%Generating Io_alpha
b_limit = 125;
Io_alpha = 1;
for k = 1:b_limit
 value_k = ((1/factorial(k))*(alpha/2).^k).^2;
 Io_alpha = Io_alpha + value_k;
end
%Generating Io_beta
Io_beta = 1;
for m = 1:b_limit
 value_m = ((1/factorial(m))*(beta/2).^m).^2;
 Io_beta = Io_beta +value_m;
end
wk_nT = Io_beta/Io_alpha;
Io_alpha
Aa
Ap
dDN
alpha
figure;
stem(n,wk_nT);
xlabel('n' );
ylabel('Amplitude' );
title('Kaiser window(Time domain)' );

%Obtaining the ideal impulse stopband filter

n_L = -(N-1)/2:1: -1;
hn_L = (1./(n_L*pi)).*(sin(low_cf*n_L*st) -sin(up_cf*n_L*st));
n_R = 1:1:(N-1)/2;
hn_R = (1./(n_R*pi)).*(sin(low_cf*n_R*st) -sin(up_cf*n_R*st));
hn_0 = 1+(2/sf).*(low_cf-up_cf);
n = [n_L,0,n_R];
h_nT = [hn_L,hn_0,hn_R];
hw_nT = h_nT.*wk_nT;
% %Plotting the ideal filter
figure;
stem(n,h_nT, ' -r' );
xlabel('n' );
ylabel('Amplitude' );
title('Ideal Impulse stopband filter(Time domain)' );

%Plotting causal stopband filter

%Plotting the noncausal stopband filter
figure;
stem(n,hw_nT);
xlabel('n' );
ylabel('Amplitude' );
title('Noncausal stopband filter window(Time domain)' );
%Plotting the causal stopband filter
shift = [0:1:N-1];
figure;
stem(shift,hw_nT);
xlabel('n' );
ylabel('Amplitude' );
title('Causal Impulse Response filter(Time Domain)' );

%Magnitude response analysis

%Question 3
[Hw,f] = freqz(hw_nT); %obtaining the frequency response and
corresponding frequencies
af = f*sf/(2*pi); %Angular frequency
log_Hw = 20.*log10(abs(Hw));
figure;
plot(af,log_Hw);
xlabel('Angular frequency(rad/s)' );
ylabel('Magnitude(dB)' );
title(' Magnitude response of the filter(Frequency domain)' );
% fvtool(af,log_Hw)
%Question 4
%Plotting the magnitude response of the passbands
%considering the lower passband
figure;
fi = round((length(af)/(sf/2)*low_cf));
hp_l = log_Hw(1:fi);
wp_l = af(1:fi);
plot(wp_l,hp_l);
axis([ -inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)' );
ylabel('Magnitude (dB)' );
title('Magnitude response of Lower Passband - Frequency Domain' );
%Considering the upperpassband
figure;
begin = round(length(af)/(sf/2)*up_cf);
wp_h = af(begin:length(af));
hp_h = log_Hw(begin:length(af));
plot(wp_h,hp_h);
axis([ -inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)' );
ylabel('Magnitude (dB)' );
title('Magnitude response of the Upper Passband - Frequency Domain' );

%Generating the input signal

%Component frequencies of the input
f1 = low_cf/2;
f2 = low_cf + (up_cf-low_cf)/2;
f3 = up_cf + (sf/2-up_cf)/2;
samples=500;
%Generating the discrete signal
n1 = 0:1:samples;
X_ax = cos(f1.*n1.*st)+cos(f2.*n1.*st)+cos(f3.*n1.*st);
figure;
subplot(2,1,1);
stem(n1,X_ax);
xlabel('n' );
ylabel('Amplitude' );
title('Input signal(Time domain)' )
subplot(2,1,2);
l_fft = 2^nextpow2(numel(n1)) -1;
x_fft = fft(X_ax,l_fft);
x_fft_plot =
[abs([x_fft(l_fft/2+1:l_fft)]),abs(x_fft(1)),abs(x_fft(2:l_fft/2+1))];
f = sf*linspace(0,1,l_fft) -sf/2;
plot(f,x_fft_plot);
xlabel('Frequency rad/s' );
ylabel('Magnitude' );
title('Input signal in the frequency domain' );
axis tight;

%Generating the output signal from our designed filter & ideal filter

%Question 6
% Filtering using frequency domain multiplication
l_fft = length(X_ax)+length(hw_nT) -1; % length for fft in x dimension
x_fft = fft(X_ax,l_fft);
hw_nT_fft = fft(hw_nT,l_fft);
out_fft = hw_nT_fft.*x_fft;
out = ifft(out_fft,l_fft);
design_out = out(floor(N/2)+1:length(out) -floor(N/2));
% Ideal Output Signal
ideal_out = cos(f1.*n1.*st)+cos(f3. *n1.*st);
%O_2 is left out because it is in the stopband
%Obtaining the output waveforms
% Frequency domain representation of output signal after filtering
using
% the designed filter
figure;
subplot(2,1,1);
l_fft = 2^nextpow2(numel(n1)) -1;
xfft_out = fft(design_out,l_fft);
x_fft_out_plot =
[abs([xfft_out(l_fft/2+1:l_fft)]),abs(xfft_out(1)),abs(xfft_out(2:l_ff
t/2+1))];
f = sf*linspace(0,1,l_fft) -sf/2;
plot(f,x_fft_out_plot);
xlabel('Frequency rad/s' );
ylabel('Magnitude' );
title('Output signal of the designed filter in the frequency domain' );
% Time domain representation of output signal after filtering using
the
% designed filter
subplot(2,1,2);
stem(n1,design_out);
xlabel('n' );
ylabel('Amplitude' );
title('Output signal of the designed filter in the time domain' );
%Obtaining the outputs of the ideal filter
figure;
subplot(2,1,1);
xfft_ideal_out = fft(ideal_out,l_fft);
x_fft_ideal_out_plot =
[abs([xfft_ideal_out(l_fft/2+1:l_fft)]),abs(xfft_ideal_out(1)),abs(xff
t_ideal_out(2:l_fft/2+1))];
plot(f,x_fft_ideal_out_plot);
xlabel('Frequency rad/s' );
ylabel('Magnitude' );
title('Output signal of the ideal filter in the frequency domain' );
% Time domain representation of output signal after filtering using
ideal filter
subplot(2,1,2);

stem(n1,ideal_out);
xlabel('n' );
ylabel('Amplitude' );
title('Output signal of the ideal filter in the time domain' );

%Comparision

mean( (design_out(:) -ideal_out(:)).^2)
mean( (design_out(:) -ideal_out(:)).^2)/(mean( (design_out(:).^2 )))
corrcoef(design_out,ideal_