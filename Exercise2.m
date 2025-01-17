%% 1. Random variables
%1.1 ROLLING A FAIR 6-FACED DICE (DISCRETE VARIABLE)
N_faces = 6; %Number of faces in the dice 
N_trials = 100000; %Number of trials (how many times the dice is rolled) 
trials = randi(N_faces,1,N_trials); % Getting random integers 1...6 


x_histogram_centers = 1:6; %x-axis for the bin center coordinates 
% x-axis for the bin edge coordinates (last edge is for the right edge) 
histogram_edges = 0.5:1:6.5; 
% Below the ~ sign means that we discard that specific output argument 
[hist_count, ~]= histcounts(trials,histogram_edges); 
figure 
bar(x_histogram_centers,hist_count) 
xlabel('Rolled number') 
ylabel('occurrence number') 
grid on % Set a grid in the plotted figure 
title('Dice rolling histogram') 
% Histograms can plotted directly as "h = histogram(trials,'BinLimits',[0.5 
% 6.5],'BinMethod','integers')", where h is the histogram object.
pdf_experimental = hist_count/N_trials;
figure 
bar(x_histogram_centers,pdf_experimental) 
xlabel('Rolled number') 
ylabel('pdf') 
grid on 
title('Dice rolling normalized histogram') 
one_face_probability = 1/N_faces; % probability of one face of the dice 
hold on % avoids removing the previous plot 
% plotting true pdf (a line between two points) 
plot([0.5 6.5],[one_face_probability one_face_probability],'r') 
hold off 
% Using legend we can name different data curves in the plot (in order of 
% appearance) 
legend('Experimental pdf','True pdf')


%1.2 NORMAL/GAUSSIAN DISTRIBUTED RANDOM VARIABLE

N_samples = 100000; 
mu = 3; 
sigma = sqrt(4); % the variance is 4 (i.e. sigma^2=4)
samples = mu+ sigma.* randn([1 N_samples]); %mean 3 and standard deviation 4

m = mean(samples);
standard_deviation = std(samples);
variance = var(samples);

bin_width = 0.5; %bin width (in the x-axis) 
bin_centers = -7:bin_width:13; %x-axis for the bin center coordinates 
% x-axis for the bin edge coordinates (last edge is for the right edge): 
% (Three dots “…” continues the command in the following line) 
bin_edges = (bin_centers(1)-bin_width/2):bin_width:(bin_centers(end)+bin_width/2); 
% ~ means that we discard that output argument 
[hist_count, ~]= histcounts(samples,bin_edges); 
pdf_experimental = hist_count/sum(hist_count*bin_width); 
figure
bar(bin_centers,pdf_experimental,1) 

xlabel('Number of samples')
ylabel('pdf')
title('Normalized histogram')

pdf_true = 1/(sqrt(2*pi)*sigma)*exp(-(mu-bin_edges).^2/(2*sigma^2)); 
hold on 
plot(bin_edges,pdf_true,'r','LineWidth',3) 

legend('Experimental pdf', 'True pdf')

b = 5.25; 
indices_with_bin_center_larger_than_b = bin_centers > b; 
considered_bin_values = pdf_experimental(indices_with_bin_center_larger_than_b); 
%area of the considered bins 
probability_X_larger_than_b = sum(considered_bin_values*bin_width);

analytic_probability = qfunc((b-mu)/sigma) ;
%% 2. Random process
%2.1 White noise vs. colored noise
N = 10000; % Number of generated samples 
noise_var = 3; % Desired noise variance 
noise = sqrt(noise_var)*randn(1,N); % Noise signal generation 

figure 
plot(noise) 
xlabel('sample index') 
ylabel('noise amplitude') 
title('White noise') 
xlim([0 100]) %define the x-axis limits 

figure 
histogram(noise,40) 
xlabel('noise amplitude') 
ylabel('histogram count') 
title('White noise histogram')
%Low-pass FIR filter
N_filter = 60; %even number 
h = firpm(N_filter,[0 0.1 0.2 1],[1 1 0 0]);
N_freq =  length(noise); 
freq_vec_filter = -1:2/N_freq:(1-2/N_freq);

figure
impz(h) %impulse response
%frequency vector values normalized between -1 and 1 
figure
plot(freq_vec_filter,10*log10(fftshift(abs(fft(h,N_freq))))) 
xlabel('Normalized frequency (F_s/2=1)') 
ylabel('Amplitude') 
title('Amplitude response of the filter')

% Filter the noise signal: 
filtered_noise = filter(h,1,noise); 
filtered_noise = filtered_noise(N_filter/2+1:end); %remove the delay 
figure 
plot(noise(1:100)) 
hold on 
plot(filtered_noise(1:100),'r') 
legend('White noise','Colored noise') 
xlabel('sample index') 
ylabel('noise amplitude') 
title('White noise and filtered (colored) noise')
figure
histogram(filtered_noise)
xlabel('Filtered noise amplitude') 
ylabel('histogram count') 
title('Filtered noise histogram')
[corr_fun, lags] = xcorr(noise); 
% we normalize the max-value to 1 and use stem-function in order to emphasize 
% the impulse-like nature of the outcome 
figure
stem(lags,corr_fun/max(corr_fun))
xlabel('\tau') 
ylabel('R(\tau)') 
title('Autocorrelation of white noise') 
xlim([-30 30])
[filtered_corr_fun, filtered_lags] = xcorr(filtered_noise); 
figure
plot(filtered_lags,filtered_corr_fun/max(corr_fun))
xlabel('\tau') 
ylabel('R(\tau)') 
title('Autocorrelation of filtered noise') 
xlim([-30 30])
noise_abs_spec = 20*log10(abs(fft(noise(1:length(filtered_noise))))); 
filtered_noise_abs_spec = 20*log10(abs(fft(filtered_noise))); 
%Define the frequency vector values (normalized between -1 and 1): 
freq_vec = -1:2/length(noise_abs_spec):1-2/length(noise_abs_spec); 
figure 
plot(freq_vec,fftshift(noise_abs_spec)) 
hold on 
plot(freq_vec,fftshift(filtered_noise_abs_spec),'r') 
hold off 
xlabel('Normalized frequency (F_s/2=1)') 
ylabel('power [dB]') 
title('Noise spectra') 
legend('White noise','Filtered (coloured) noise') 

%2.2 Random walk model

N_samples = 2000; %Number of samples for each realization 
N_ensemble = 5000; %Number of signal realizations (i.e., the size of ensemble) 
%Step probability and step size: 
p = 0.5; % P(Wi=s) = p, P(Wi=-s) = 1-p 
s = 1; %step length 
n = 1:N_samples; % vector of samples indices 
% Generating matrix of randomly generated steps: 
W = rand(N_ensemble,N_samples); 
% (i.e. uniformly distributed random values between 0 and 1) 
indices_with_positive_s = W<p; % find out steps going "up" 
W(indices_with_positive_s) = s; % Define steps for going "up" 
W(~indices_with_positive_s) = -s; % Define steps for going "down" 

% The overall "random walk" is achieved by taking the cumulative sum over the 
% steps: 
X = cumsum(W,2);  

figure 
for ind = 1:5 
    subplot(5,1,ind) 
    plot(n,X(ind,:)) 
    xlabel('n') 
    ylabel('X(n)') 
    grid on 
    title(['Realization #' num2str(ind)]) 
   
end 
 
set(gcf,'units','normalized','outerposition',[0 0 1 1]) 
mean_theory = n*s*(2*p-1); % Theoretical mean 
var_theory = n*(2*s)^2*p*(1-p); % Theoretical variance 
mean_observed = mean(X); % Empirical mean  
var_observed = var(X); % Empirical variance 
figure 
plot(n,mean_observed,'b','LineWidth',3) 
hold on 
plot(n,mean_theory,'r:','LineWidth',2) 
hold off 
legend('observed mean','theoretical mean') 
ylim([-2 2]) % set the axis limits in y-direction only 
xlabel('n') 
ylabel('Mean') 
title('Mean over the sample index') 


figure
plot(n,var_observed,'b','LineWidth',3)
hold on 
plot(n,var_theory,'r:','LineWidth',2) 
hold off 
legend('observed variance','theoretical variance') 

ylim([0 20]) 

xlabel('n') 
ylabel('Variance') 
title('Variance over the sample index') 
