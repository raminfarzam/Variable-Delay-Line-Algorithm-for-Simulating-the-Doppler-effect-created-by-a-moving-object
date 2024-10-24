%% The program is primirally is desigend for eimulating the effect of the moving object and
clc; clear all; close all
% Audio signal from that should be played from one roatating sound Sources
par.audioFilePath = 'Miley Cyrus - Wrecking Ball Vocals Only (320 kbps)(1) (mp3cut.net).flac';
[par.audioData, par.sampleRate] = audioread(par.audioFilePath);
% The selected audio has same left and right channel signal so there is no
% relative delay
par.leftChannel  = par.audioData(:, 1);
par.rightChannel = par.audioData(:, 1);
%% Display original audio
figure(11);
subplot(1, 1, 1);
plot((1:length(par.leftChannel))/par.sampleRate, par.leftChannel, 'r*');
title('Left Channel');hold on
xlim([0 , length(par.leftChannel)/par.sampleRate])
subplot(1, 1, 1);
plot((1:length(par.rightChannel))/par.sampleRate, par.rightChannel, 'b');
title('Right Channel');hold off
%% Down Sample the Signal with the Factor of 2,4,8  
fs = par.sampleRate./ 2 ;           % down sampeling frequency
s1 = resample(par.leftChannel,fs , par.sampleRate);
s2 = resample(par.rightChannel,fs , par.sampleRate);
Ns = length(s1);
%% generate the sound at each of the Microphone
par.d = 0.18 ;                              % location of the Microphones (two ears of a person)
par.omega = pi/10;                          % w = pi [rad/s]
par.r0 = 5 ;                              	% 3[m]
par.N  = 2;
par.v = 335;                                        % speed of wavefront 
t = (0:Ns-1)./fs;                                   % time axis
par.x = par.r0 .* cos(par.omega * t);               % x speaker
par.y = par.r0 .* sin(par.omega * t);               % y speaker
theta_true = linspace (0 ,2*pi , Ns);               % tetha speaker
r1  =  sqrt( par.y.^2 + (par.x + par.d/2).^2 ) ;	% r1
r2  =  sqrt( par.y.^2 + (par.x - par.d/2).^2 ) ;	% r2
DTOA_true = (cos( theta_true ).* par.d ./ par.v)';          
DTOA_true1 = ( r1 - r2 )./par.v;
delayl = (r1./ par.v ).*fs;                         % delay observed at left microphone
delayr = (r2./ par.v ).*fs;                         % delay observed at right microphone
% remove the base deay that is due to distance and onl consider the doppler
% effect delays
d_sample_l = round( delayl - min(delayl) );
d_sample_r = round( delayr - min(delayr) );
plot(d_sample_l);
ylabel('delay sample'); xlabel('sample number')
title('variable delay due to doppler effect of the circular motion ')
%% variable dealy line Algorithm
s1_rot = Variable_Delay(s1,d_sample_l);     %use the delay line function and only consider the linear intelpolation
s2_rot = Variable_Delay(s2,d_sample_r);
s1_rot =  (par.r0 ./ r1).^2 .* s1_rot ;     % attenuated signal at the right mic
s2_rot =  (par.r0 ./ r2).^2 .* s2_rot ;     % attenuated signal at the left mic
s1_rot = s1_rot'; s2_rot = s2_rot';         
par.y = [s1_rot , s2_rot];
%% Play the rotated sound
figure(22)
plot(t,s1_rot , 'r');hold on
plot(t,s2_rot , 'b');
xlim([0 t(end)])
xlabel('Time (s)'); ylabel('Amplitude')
title('Recieved signal at each michrophone')
% Define the region to magnify
magnifyRegionX = [1.10, 1.135];
magnifyRegionY = [-0.55, 0.7];
% Create a dashed line box around the magnified region
rectangle('Position', [magnifyRegionX(1), magnifyRegionY(1), ...
                       diff(magnifyRegionX), diff(magnifyRegionY)], ...
          'EdgeColor', 'green', 'LineStyle', '--');
% Create an arrow to indicate the magnified region
arrowStart = [0.17, 0.4];
arrowEnd = [0.3, 0.75];
annotation('arrow', arrowStart, arrowEnd, 'Color', 'green');
% Create a magnified subplot near the arrow
magnifiedPlot = axes('Position', [0.2, 0.75, 0.2, 0.1]);
box on;
plot(t,s1_rot, 'b-');hold on
plot(t,s2_rot,'r-');hold on
xlim(magnifyRegionX);
ylim(magnifyRegionY);
title('Magnified Plot');
magnifyRegionX = [12.10, 12.135];
magnifyRegionY = [-0.55, 0.7];
% Create a dashed line box around the magnified region
rectangle('Position', [magnifyRegionX(1), magnifyRegionY(1), ...
                       diff(magnifyRegionX), diff(magnifyRegionY)], ...
          'EdgeColor', 'green', 'LineStyle', '--');
% Create an arrow to indicate the magnified region
arrowStart = [0.17, 0.4];
arrowEnd = [0.3, 0.75];
annotation('arrow', arrowStart, arrowEnd, 'Color', 'green');
% Create a magnified subplot near the arrow
magnifiedPlot = axes('Position', [0.5, 0.75, 0.3, 0.1]);
box on;
plot(t,s1_rot, 'b-');hold on
plot(t,s2_rot,'r-');hold on
xlim(magnifyRegionX);
ylim(magnifyRegionY);
title('Magnified Plot');
% sound([s1_rot , s2_rot],fs);
% player = audioplayer(par.y,fs);
% play(player);
%% Save the Audio for the output
% filename = 'Audio_for_pitch.wav';  
% audiowrite(filename, par.y, fs);
%% Add Noise
Noise.SNR_dB = -5:5:30;                         % SNR levels in dB
Noise.numSNRLevels = length(Noise.SNR_dB);      % SNR Length
y_noisy = zeros(Ns , 2.*Noise.numSNRLevels);
for n = 1: Noise.numSNRLevels                   % add noise at different SNR levels
    snr = Noise.SNR_dB(n);
    noisePower = 10^(- snr / 10); 
    noise1 = noisePower * randn(Ns, 1);
    noise2 = noisePower * randn(Ns, 1);
    % Add noise to the original signal
    noisy1 = par.y(:,1) + noise2;
    noisy2 = par.y(:,2) + noise1;
    y_noisy(:, 2*n-1) = noisy1;
    y_noisy(:, 2*n) = noisy2;
end
%%  Find Maximum corelation for the different intervals that &τ(t) can be considered constant
Windows =[50, 100, 200 ,500 ,1000 ,1500, 5000];           % different # of windows
DTOA_est = zeros( Ns , length(Windows) );
MSE = zeros(1, length(Windows));
for w = 1:length(Windows)               % loop over different windows length
    w_len = Windows(w);
    N = round( Ns / w_len);
    tau_est = zeros(w_len, 1);
    for k = 1:w_len                         % loop over w_len 
%         sliced1 = s1_rot((k - 1) * N + 1 : k * N);
%         sliced2 = s2_rot((k - 1) * N + 1 : k * N);
        sliced1 = y_noisy((k - 1) * N + 1 : k * N,end-1);
        sliced2 = y_noisy((k - 1) * N + 1 : k * N,end);
        convolution = conv(sliced1, sliced2(end : -1:1)); 
        [~, mind] = max(convolution);
        
        if mind == 1 || mind == 2*N - 1
            tau_est(k) = NaN;
        else              % quadratic interpolation
            delta = 0.5 * (convolution(mind - 1) - convolution(mind + 1)) / (convolution(mind - 1) + convolution(mind + 1) - 2 * convolution(mind));
            tau_est(k) = (delta + mind - N) / fs;
        end
    end

    DTOA_est(:, w) = interp1(1:w_len, tau_est, linspace(1, w_len, Ns));
    err = DTOA_est(:, w) - DTOA_true;
    err(isnan(err)) = [];
    MSE(w) = err' * err / w_len;
end
[~, best] = min(MSE);                                                      % choose the best block size
%% iteration for different SNRs with proper choice of the Window
DTOA_est_SNR = zeros( Ns , Noise.numSNRLevels );
MSE_SNR = zeros(1, Noise.numSNRLevels);
w_len = Windows(best);
N = round( Ns / w_len);
for  i =1:Noise.numSNRLevels
    tau_est = zeros(w_len, 1);
    for k = 1:w_len                         % loop over w_len 
        sliced1 = y_noisy((k - 1) * N + 1 : k * N, (2*i) - 1);
        sliced2 = y_noisy((k - 1) * N + 1 : k * N, 2*i );
        convolution = conv(sliced1, sliced2(end : -1:1)); 
        [~, mind] = max(convolution);
        % quadratic interpolation
        if mind == 1 || mind == 2*N - 1
            tau_est(k) = NaN;
        else
            delta = 0.5 * (convolution(mind - 1) - convolution(mind + 1)) / (convolution(mind - 1) + convolution(mind + 1) - 2 * convolution(mind));
            tau_est(k) = (delta + mind - N) / fs;
        end
    end
    DTOA_est_SNR(:, i) = interp1(1:w_len, tau_est, linspace(1, w_len, Ns));
    err = DTOA_est_SNR(:, i) - DTOA_true;
    err(isnan(err)) = [];
    MSE_SNR(i) = err' * err / w_len;
    CRB(i) = (3* (1 + 2*10^( Noise.SNR_dB(i) / 10)))/(pi^2 * N *8* 10^( 2*Noise.SNR_dB(i) / 10));
end
%%
figure(33);
subplot(2,1,1)
plot(Windows, MSE ,'*-');hold on
title(['MSE vs different Number of windows in SNR =' ,num2str(Noise.SNR_dB(best)), 'dB']);
xlabel(' Number of Windows');
ylabel('MSE (Mean Squared Error)');grid on;

% figure(44);
subplot(2,1,2)
plot(t,DTOA_true);hold on;plot(t, DTOA_est(:, best))
xlabel('time(s)'); ylabel('\Delta\tau (t)'); xlim([0 t(end)]);
title(['Estimation of the \Delta \tau(t) for the Best window size, and SNR =' ,num2str(Noise.SNR_dB(best)), 'dB'])
legend('true value', 'estimated \Delta \tau(t) for the best number of windows')
grid on;
theta_est = acos(par.v / (par.d) * DTOA_est(:, best));              % from DToA to DoA
theta_TRUE = acos(par.v / (par.d) * DTOA_true);
figure(55);
plot(abs(theta_est));hold on;
plot(theta_TRUE)
legend('true value', 'estimated \Delta \tau(t) for the best number of windows')
grid on;
%%
figure(66);
semilogy(Noise.SNR_dB, MSE_SNR ,'*-');hold on; plot(Noise.SNR_dB, CRB,'r^-')
title(['MSE of \Delta \tau(t) estimation vs different SNR for # of windows = ',num2str(Windows(best)), 's']);
xlabel(' SNR (dB)');
ylabel('MSE (Mean Squared Error)');grid on;
legend('MSE', 'CRB')

% figure(77);
% plot(t,DTOA_true);hold on;
% plot(t, DTOA_est_SNR(:,3:end));
% xlabel('time(s)'); ylabel('\Delta\tau (t)'); xlim([0 t(end)]);ylim([min(DTOA_true) max(DTOA_true)])
title('Estimation of the \Delta \tau(t) for different SNR')
% legend('true value', 'estimated \Delta \tau(t) for the best window intervals')
grid on;
% theta_est = acos(par.v / (par.d) * DTOA_est(:, best)); % from DToA to DoA
% figure(55);
% plot(abs(theta_est));hold on;
% plot(theta_true)
% grid on;
%% part B realigh the Signals 
a_1 = zeros(1, Ns);
a_2 = zeros(1, Ns);
do1 = DTOA_est(:, best).*fs; 
a_1(do1 < 0) = -do1(do1 < 0) ;
a_2(do1 > 0) =  do1(do1 > 0)  ;
a_1 = floor(real(a_1));
a_2 = floor(real(a_2));
% figure();plot(a_1);
% hold on; plot(a_2)
%% REconstruct R_1 and R_2 
R_1 = Variable_Delay(y_noisy(:,end-1),a_1);
R_2 = Variable_Delay(y_noisy(:,end),a_2);
%%
% plot(R_1);hold on;
% plot(R_2);
hold on;
pp =[zeros(max(a_1),1);s1];
plot(pp);hold on
legend('R_1' , 'R_2', 'original shifted')
% Define the region to magnify
magnifyRegionX = [142600, 142785];
magnifyRegionY = [-0.55, 0.7];
% Create a dashed line box around the magnified region
rectangle('Position', [magnifyRegionX(1), magnifyRegionY(1), ...
                       diff(magnifyRegionX), diff(magnifyRegionY)], ...
          'EdgeColor', 'green', 'LineStyle', '--');
% Create an arrow to indicate the magnified region
% arrowStart = [0.17, 0.4];
% arrowEnd = [0.3, 0.75];
% annotation('arrow', arrowStart, arrowEnd, 'Color', 'green');
% Create a magnified subplot near the arrow
magnifiedPlot = axes('Position', [0.2, 0.75, 0.2, 0.1]);
box on;
plot(R_1, 'b-');hold on
plot(R_2,'r-');hold on
plot(pp,'y-');hold on
xlim(magnifyRegionX);
ylim(magnifyRegionY);
title('Magnified Plot');
% player = audioplayer([R_1 , R_2],fs);
% play(player);
R_hat =  0.5.*R_1 + 0.5.*R_2;
mse_result = calculate_mse(R_hat, pp');
player = audioplayer([R_hat' , pp(1:Ns)],fs);
% play(player);
disp(['Mean Squared Error: ', num2str(mse_result)]);
%% part C realignement for the different SNRs
a_1 = zeros(Noise.numSNRLevels, Ns);
a_2 = zeros(Noise.numSNRLevels, Ns);
for nn = 1: Noise.numSNRLevels
    do1 = DTOA_est_SNR(:, nn).*fs; 
    a_1(nn,do1 < 0) = -do1(do1 < 0) ;
    a_2(nn,do1 > 0) =  do1(do1 > 0)  ;
    a_1(nn,:) = floor(real(a_1(nn,:)));
    a_2(nn,:) = floor(real(a_2(nn,:)));
    % REconstruct R_1 and R_2 
    R_1 = Variable_Delay( y_noisy(:,2* nn -1) , a_1(nn,:));
    R_2 = Variable_Delay( y_noisy(:,2* nn ) , a_2(nn,:));

    pp =[zeros(max(a_1(nn,:)),1);s1];
    
    R_hat =  0.5.*R_1 + 0.5.*R_2;
    mse_result_SNR(nn) = calculate_mse(R_hat, pp');
    RR =[zeros(min(d_sample_l) + floor(mean(a_1(nn,:))),1);s1]';
    tau_est = zeros(w_len, 1);
    
    for k = 1:w_len                         % loop over w_len 
        sliced1 = R_hat((k - 1) * N + 1 : k * N);
        sliced2 = RR((k - 1) * N + 1 : k * N);
        convolution = conv(sliced1, sliced2(end : -1:1)); 
        [~, mind] = max(convolution);
        % quadratic interpolation
        if mind == 1 || mind == 2*N - 1
            tau_est(k) = NaN;
        else
            delta = 0.5 * (convolution(mind - 1) - convolution(mind + 1)) / (convolution(mind - 1) + convolution(mind + 1) - 2 * convolution(mind));
            tau_est(k) = (delta + mind - N) / fs;
        end
    end
    ToA_est(nn,:) = interp1(1:w_len, tau_est, 1:w_len);
    err = ToA_est(nn,:) - (par.r0 / par.v).*ones(1,w_len);
    err(isnan(err)) = [];
    MSE_TOA(nn) = err * err' / w_len;
    
%     correlation = xcorr(R_hat' , RR);
%     % Find the time shift
%     [~, idx] = max(correlation);
%     time_shift = idx - Ns + 1;
%     time_of_arrival = time_shift / fs;
%     MSE_TOA(nn) = time_of_arrival - (par.r0 / par.v);
end
%%
figure('Name','MSE_part_c')
semilogy(Noise.SNR_dB, mse_result_SNR,'-*'); grid on
title('MSE of $$\hat{R(t)}$$ vs the Original Signal','Interpreter','latex')
xlabel('SNR(dB)')
ylabel('MSE')
figure('Name','TOA')
plot(Noise.SNR_dB, MSE_TOA) 
title('MSE of TOA estimatoion Using $$\hat{R(t)}$$','Interpreter','latex')
xlabel('SNR(dB)')
ylabel('MSE');grid on
