% BASIC OFDM CODE
% UNCODED QPSK - OVERALL RATE IS 2
close all
clear all
clc
SNR_dB = 40;% SNR PER BIT
NUM_FRAMES = 10^2; 
FFT_LEN = 1024;
NUM_BIT = 2*FFT_LEN; % NUMBER OF DATA BITS
CHAN_LEN = 10; % NUMBER OF CHANNEL TAPS
CP_LEN = CHAN_LEN-1; % LENGTH OF THE CYCLIC PREFIX
FADE_VAR_1D = 0.5; % 1D FADE VARIANCE OF THE CHANNEL
FADE_STD_DEV = sqrt(FADE_VAR_1D); % STANDARD DEVIATION OF THE FADING CHANNEL
% SNR PER BIT DEFINITION - OVERALL RATE IS 2
SNR = 10^(0.1*SNR_dB); % LINEAR SCALE
NOISE_VAR_1D = 0.5*2*2*CHAN_LEN*FADE_VAR_1D/(2*SNR*FFT_LEN); % 1D AWGN VARIANCE 
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION
tic()
C_BER = 0; % bit errors in each frame
for FRAME_CNT = 1:NUM_FRAMES
%----            TRANSMITTER      -----------------------------------------
% SOURCE
A = randi([0 1],1,NUM_BIT);
% QPSK MAPPING
F_SIG_NO_CP = 1-2*A(1:2:end) + 1i*(1-2*A(2:2:end));
% IFFT 
T_SIG_NO_CP = ifft(F_SIG_NO_CP);
% INSERTING CYCLIC PREFIX
T_SIG_CP = [T_SIG_NO_CP(end-CP_LEN+1:end) T_SIG_NO_CP];
%---------------     CHANNEL      -----------------------------------------
% RAYLEIGH FREQUENCY SELECTIVE FADING CHANNEL
FADE_CHAN = normrnd(0,FADE_STD_DEV,1,CHAN_LEN)+1i*normrnd(0,FADE_STD_DEV,1,CHAN_LEN);
FREQ_RESP = fft(FADE_CHAN,FFT_LEN); % ACTUAL CHANNEL FREQUENCY RESPONSE
% AWGN
AWGN = normrnd(0,NOISE_STD_DEV,1,FFT_LEN+CP_LEN+CHAN_LEN-1)+1i*normrnd(0,NOISE_STD_DEV,1,FFT_LEN+CP_LEN+CHAN_LEN-1);
% CHANNEL OUTPUT
T_REC_SIG = conv(T_SIG_CP,FADE_CHAN) + AWGN;
%----------------      RECEIVER  ------------------------------------------
% CP & TRANSIENT SAMPLES REMOVAL
T_REC_SIG(1:CP_LEN) = [];
T_REC_SIG_NO_CP = T_REC_SIG(1:FFT_LEN);
% PERFORMING THE FFT
F_REC_SIG_NO_CP = fft(T_REC_SIG_NO_CP);
% ML DETECTION
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
QPSK_SYM1 = QPSK_SYM(1)*ones(1,FFT_LEN);
QPSK_SYM2 = QPSK_SYM(2)*ones(1,FFT_LEN);
QPSK_SYM3 = QPSK_SYM(3)*ones(1,FFT_LEN);
QPSK_SYM4 = QPSK_SYM(4)*ones(1,FFT_LEN);
DIST = zeros(4,FFT_LEN);
DIST(1,:)=abs(F_REC_SIG_NO_CP - FREQ_RESP.*QPSK_SYM1).^2; 
DIST(2,:)=abs(F_REC_SIG_NO_CP - FREQ_RESP.*QPSK_SYM2).^2;
DIST(3,:)=abs(F_REC_SIG_NO_CP - FREQ_RESP.*QPSK_SYM3).^2;
DIST(4,:)=abs(F_REC_SIG_NO_CP - FREQ_RESP.*QPSK_SYM4).^2; 
% COMPARING EUCLIDEAN DISTANCE
[~,INDICES] = min(DIST,[],1);
% MAPPING INDICES TO QPSK SYMBOLS
DEC_QPSK_MAP_SYM = QPSK_SYM(INDICES);
% DEMAPPING QPSK SYMBOLS TO BITS
DEC_A = zeros(1,NUM_BIT);
DEC_A(1:2:end) = real(DEC_QPSK_MAP_SYM)<0;
DEC_A(2:2:end) = imag(DEC_QPSK_MAP_SYM)<0;
% CALCULATING BIT ERRORS IN EACH FRAME
C_BER = C_BER + nnz(A-DEC_A);
end
toc()
BER = C_BER/(NUM_BIT*NUM_FRAMES)