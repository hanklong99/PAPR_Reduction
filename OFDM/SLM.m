clear;
close all;
clc

%% Initialization
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

%Original
bers = zeros(size(EbN0_dB));
Power = zeros(Nsym*Nframe,N_iter);
Average=zeros(N_iter,1);
Peak=zeros(N_iter,1);
%SLM
N_SLM = 20; % number of sets of SLM Phase
bers_SLM=zeros(size(EbN0_dB));
Power_SLM=zeros(N_SLM, Nsym*Nframe);
OFDM_SLM=zeros(N_SLM,Nsym*Nframe);
PAPR_SLM=zeros(N_SLM,N_iter);
min_PAPR = zeros(N_iter,1);


%% 1??Calculate the amplitude of Noise
% Alternative method
EbN0=10.^(EbN0_dB/10);
%1/64 is the symbol energy
sigPow = 1/64;
Eb=sigPow/Nbps*T;  %bit energy is half of symbol energy
N0=Eb./EbN0; %sqrt(2) for normalizing the noise magnitude


%% 2??Generate signal again
N_block = N_iter;
%gain = 1/sqrt(2)*[randn(1, N_block) + 1j*randn(1, N_block)];
for k_EbN0 = 1:length(EbN0_dB)
    
    ber_count = 0; % error bit number
    ber_count_SLM=0;
    for m = 1:N_iter
        
        symbol = randi(N_LEVEL,symbol_total,1) - 1; % informatio source
        % Info Source & Mod
        graycode = symbol;
        % Gray code (according to website)
        index_3 = find(symbol==3);
        index_2 = find(symbol==2);
        graycode(index_3) = 2;
        graycode(index_2) = 3;
        % mapped QPSK data
        data_QPSK = exp(1j * graycode* pi/2 + 1j * pi/4);
        ofdm=OFDM_Gen(data_QPSK,Ng, Nframe, Nfft, Nsym);
        Power(:,m)=abs(ofdm).^2;
        % total energy        %Power(:,m)=abs(x_GI).^2;
        Average(m)=mean(Power(:,m));
        Peak(m)=max(Power(:,m));


        % Channel signal
        y = ofdm;
        %AWGN noise
        % sqrt(2) normalization
        complex_noise = sqrt(N0(k_EbN0))*(randn(size(y))+1j*randn(size(y)))/sqrt(2); 
        %Add Noise
        Rx=y+complex_noise;
        %Calculate BER
        SLM_Phase=zeros(length(data_QPSK), 1); % Clipping do not need SLM Phase
        ber_count_temp = OFDM_Con(SLM_Phase,Rx, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
        ber_count=ber_count+ber_count_temp;

         % Selective Mapping
         gray_P = randi([0 N_LEVEL-1],length(data_QPSK), N_SLM);
         gray_P=exp(1j * gray_P* pi/2 + 1j * pi/4);
         for i=1:N_SLM
             QPSK_SLM = data_QPSK.*gray_P(:,i);
             OFDM_SLM(i,:) = OFDM_Gen(QPSK_SLM,Ng, Nframe, Nfft, Nsym);
             Power_SLM(i,:) = abs(OFDM_SLM(i,:)).^2;
             PAPR_SLM(i,m)=(max(Power_SLM(i,:))/mean(Power_SLM(i,:)));
         end
         min_PAPR(m)=min(PAPR_SLM(:,m));
         [a,b]=find(PAPR_SLM==min_PAPR(m));
         index_min=a;
         Rx_SLM = OFDM_SLM(index_min, :) + complex_noise;
         ber_temp_SLM = OFDM_Con(gray_P(:,index_min),Rx_SLM, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
         ber_count_SLM = ber_count_SLM+ber_temp_SLM;

    end
    counter=[k_EbN0 length(EbN0_dB)];
    disp(counter);
    % 2 bits per symbol
    ber = ber_count/symbol_total/N_iter/2;
    bers(k_EbN0) = ber;

    ber_SLM=ber_count_SLM/symbol_total/N_iter/2;
    bers_SLM(k_EbN0)= ber_SLM;
end

%% 4??Plot
M = N_LEVEL;
Pe_AWGN = ber_QAM(0:EbN0_dB(end),M,'AWGN');
%Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'Rayleigh');
figure(2);
% AWGN Analytic
semilogy(0:EbN0_dB(end),Pe_AWGN,'r-');hold on;grid on;
semilogy(EbN0_dB,bers,'b-*');% Simulation
semilogy(0:EbN0_dB(end),bers_SLM, 'p-');% SLM Simulation

legend('AWGN Analytic', 'AWGN Simulation', 'SLM Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 15 1e-7 1]);

figure(1);

papr_db=10*log10(Peak./Average);
[f,x]=ecdf(papr_db);
ccdf=1-f;

papr_db_SLM=10*log10(min_PAPR);
[f_SLM,x_SLM]=ecdf(papr_db_SLM);
ccdf_SLM=1-f_SLM;

marker_idx = 1:1e3:length(x);
semilogy(x,ccdf,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_SLM,ccdf_SLM,'*-r','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
%legend('CCDF of Original', 'CCDF of Clipping','CCDF of SLM');
legend('CCDF of Original', 'CCDF of Clipping');