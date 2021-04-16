%% QPSK-OFDM

% Add Cyclic Prefix
% QPSK using Gray Code
% AWGN Channel

% 1、Calculate the power（loop for N_iter times）
% 2、Calculate the noise power(according to signal power)
% 3、Generate signal and noise agian
% 4、Decide and calculate error bits
% repeat N_iter times

clear;
close all;
clc

EbN0_dB = 0:1:15;
N_iter = 1e4; % repeat N_iter
% N_iter = 1e5;
A=1;
T=1;
Nbps = 2;
N_LEVEL = 2^Nbps; % QPSK levels
Nfft = 64;        % numbers of fft and ifft
Nvc = 0;          % virtual carrier
Nused = Nfft - Nvc; % No virtual carrier
Ng = Nfft/4;      % Guard interval
Nsym = Nfft + Ng; % total number of symbols per frame
Nframe = 3;       % the number of frames
symbol_total = Nfft * Nframe; % information source
N_SLM = 10;
%initialize the arrays
Power_SLM=zeros(N_SLM, Nsym*Nframe);
OFDM_SLM=zeros(N_SLM,Nsym*Nframe);
PAPR_SLM=zeros(N_SLM,N_iter);
min_PAPR = zeros(N_iter,1);
x_GI = complex(zeros(Nsym * Nframe,1),0);
Y = complex(zeros(symbol_total,1),0);
symbol_output = zeros(symbol_total,1);
bers = zeros(size(EbN0_dB));
bers_clip = zeros(size(EbN0_dB));
bers_SLM=zeros(size(EbN0_dB));
Power = zeros(Nsym*Nframe,N_iter);
Average=zeros(N_iter,1);
Peak=zeros(N_iter,1);
Power_clip = zeros(Nsym*Nframe,N_iter);
Average_clip=zeros(N_iter,1);
Peak_clip=zeros(N_iter,1);
mu_A = 1/8;
mu = exp(1)-1;
%% 1、Calculate the signal power
sigEng = 0; %signal energy
for m = 1:N_iter
    symbol = randi(N_LEVEL,symbol_total,1) - 1; % information source
    % Info Source & Mod
    graycode = symbol;
    % Gray code (according to website)
    index_3 = find(symbol==3);
    index_2 = find(symbol==2);
    graycode(index_3) = 2;
    graycode(index_2) = 3;
    % mapped QPSK data
    data_QPSK = exp(1j * graycode* pi/2 + 1j * pi/4);
     % Mu Companding
         OFDM_Compd = OFDM_Gen(data_QPSK,Ng, Nframe, Nfft, Nsym);   
         mu_QPSK= (mu_A*log(1+mu*abs(OFDM_Compd)/mu_A)/log(1+mu)).*sign(OFDM_Compd);
         Power(:,m)=abs(mu_QPSK).^2;
    sigEng = sigEng + sum(abs(mu_QPSK)); % conjugate to calculate energy
end

% signal power of one symbol
sigPow = sigEng/Nsym/Nframe/N_iter;
disp(sigPow);
