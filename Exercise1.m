%% Task 1
Fs=500; % sampling frequency = 10*50Hz 
 Ts=1/Fs; % sampling interval 
 t=0:Ts:0.1; % sampling time instants [s] 
 x=sin(2*pi*50*t); % signal vector 
 figure % Opens a new figure 
 plot(t,x) % plot in time domain 
 %% Task 2

F_x=fft(x); % DFT of X, saved to Fx 
Nx=length(x); 
Fo=1/(Ts*Nx); % frequency resolution. 
freq1=0:Fo:(Nx-1)*Fo; % One-sided frequency Axis 
figure 
plot(freq1,abs(F_x)/Nx) % One-sided amplitude Spectrum 
freq2=-Nx/2*Fo:Fo:(Nx/2-1)*Fo; % Two-sided frequency Axis 
figure 
plot(freq2,fftshift(abs(F_x)/Nx)) % Two-sided amplitude Spectrum
%% Task 3
clear;
%3.1
Fs = 16e9; % sampling frequency Fs = 16 GHz
Ts = 1/Fs; %sampling intervals
T=20e-9;
t=0:Ts:T;
N = length(t);
fc = 800e6; % signal's frequency
x = sin(2*pi*fc*t);

figure
plot(t,x)
title('Time domain Plot of x(t)') 
xlabel('t [s]') 
ylabel('Amplitude') 
axis([0 20e-9 -1.2 1.2])

Fo =1/(Ts*N); %frequency resolution
Fx = fft(x); %DFT
f = -N/2*Fo:Fo:Fo*(N/2-1);% two-sided frequency axis

figure
plot(f,fftshift(abs(Fx)/N))

%3.2
m=sin(2*pi*750e6*t);
s=x.*m;

figure
plot(t,s)

Fx = fft(s);
figure
plot(f,fftshift(abs(Fx)/N))

% 3.3
n=randn(size(x));
y=10*s+n;

figure
plot(t,y)

F_y = fft(y);
figure
plot(f/1e6,fftshift(abs(F_y))/N) 

%% Task 4
%4.1
f_cut = 200e6;
order = 10;
fr = f_cut/(Fs/2); % Cut-off frequency normalized to 1. 
[b,a] = butter(order,fr);% Coefficients of the filter 

freqz(b,a,N,Fs)
title('Frequency response of the Butterworth filter')

y_filtered_butter = filter(b,a,y);
figure
plot(t,y_filtered_butter)

F_y_filtered = fft(y_filtered_butter);
figure
plot(f/1e6,fftshift(abs(F_y_filtered))/N)


% 4.2
order = 60;
f_filter = [0 0.8e9 1.3e9 1.8e9 2.3e9 Fs/2]/(Fs/2); 
a_filter = [0 0 1 1 0 0]; 
b = firpm(order,f_filter,a_filter); 

figure
stem(-order/2:order/2,b) 

F_b = fft(b,N);
figure
plot(f/1e6,fftshift(abs(F_b))/N)

y_filtered_FIR = filter(b,1,y);

figure
plot(t,y_filtered_FIR)

F_y_fir = fft(y_filtered_FIR);
figure
plot(f/1e6,fftshift(abs(F_y_fir))/N)
