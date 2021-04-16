clear;
close all;
clc

%% Initialization
EbN0_dB = 0:1:30;
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
%TSLM
N_TSLM = 10; % number of sets of SLM Phase
Power_TSLM=zeros(N_TSLM, Nsym*Nframe);
PAPR_TSLM=zeros(N_TSLM, N_iter);
OFDM_TSLM=zeros(N_TSLM,Nsym*Nframe);
min_PAPR = zeros(N_iter,1);
T_TSLM = zeros(Nfft,Nfft, N_TSLM);
F = zeros(Nfft,Nfft, N_TSLM);
s_GI = zeros(Nsym*Nframe,1);
%Comapnding
bers_Compd_TSLM=zeros(size(EbN0_dB));
Power_Compd_TSLM=zeros(Nsym*Nframe, N_iter);
Average_Compd_TSLM=zeros(N_iter,1);
Peak_Compd_TSLM=zeros(N_iter,1);
OFDM_Compd=zeros(Nsym*Nframe, N_iter);

mu=5;
mu_A = 2;

%% 1 Calculate the amplitude of Noise
% Alternative method
EbN0=10.^(EbN0_dB/10);
%1/64 is the symbol energy
sigPow = 1/64;
Eb=sigPow/Nbps*T;  %bit energy is half of symbol energy
N0=Eb./EbN0; %sqrt(2) for normalizing the noise magnitude


%% 2 Generate signal again
N_block = N_iter;
%gain = 1/sqrt(2)*[randn(1, N_block) + 1j*randn(1, N_block)];
for k_EbN0 = 1:length(EbN0_dB)
    
    ber_count = 0; % error bit number
    ber_count_Compd_TSLM=0;
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
        SLM_Phase=zeros(length(data_QPSK), 1); % do not need SLM Phase
        ber_count_temp = OFDM_Con(SLM_Phase,Rx, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
        ber_count=ber_count+ber_count_temp;
        
         % Mu Companding
         OFDM_Compd = OFDM_Gen(data_QPSK,Ng, Nframe, Nfft, Nsym);   
                  % Time-domain Selective Mapping
         gray_P = randi([0 N_LEVEL-1],length(data_QPSK)/Nframe, N_TSLM);
         gray_P=exp(1j * gray_P* pi/2 + 1j * pi/4);
         %DFT Matrix
         D = dftmtx(Nfft);
         Inv_D= inv(D);
         % Phase rotation
         for i=1:N_TSLM
             F(:,:,i)=diag(gray_P(:,i));
             T_TSLM(:,:,i)=Inv_D*F(:,:,i)*D;
             %S=F.*data_QPSK;
             %s=Inv_D*S.';
             %t_QPSK=Inv_D*data_QPSK;
             % Add guard interval
             for k=1:Nframe
             % effective symbols
             kk1 = (k-1) * Nfft + (1:Nfft);
             % total symbols
             kk2 = (k-1) * Nsym + (1:Nsym);
             t_QPSK = Inv_D*data_QPSK(kk1);
             s = T_TSLM(:,:,i)*t_QPSK;
             s_GI(kk2) = guard_interval(Ng,Nfft,s);
             end
             OFDM_TSLM(i,:)=s_GI;
             Power_TSLM(i,:) = abs(OFDM_TSLM(i,:)).^2;
             PAPR_TSLM(i,m)=(max(Power_TSLM(i,:))/mean(Power_TSLM(i,:)));
         end
         
         min_PAPR(m)=min(PAPR_TSLM(:,m));
         [a,b]=find(PAPR_TSLM==min_PAPR(m));
         index_min=a;
         mu_QPSK= (mu_A*log(1+mu*abs(OFDM_TSLM(index_min, :))/mu_A)/log(1+mu)).*sign(OFDM_TSLM(index_min, :));
         Power_Compd_TSLM(:,m)=abs(mu_QPSK).^2;
         % total energy
         Average_Compd_TSLM(m)=mean(Power_Compd_TSLM(:,m));
         Peak_Compd_TSLM(m)=max(Power_Compd_TSLM(:,m));
         % The signal has been compressed, calculate new noise
         complex_noise_Compd = sqrt(N0(k_EbN0)*64*0.0445)*(randn(size(y))+1j*randn(size(y)))/sqrt(2); 
         %Received
         Rx_Compd_TSLM = mu_QPSK + complex_noise_Compd;
         %Decompanding
         Rx_DeCompd_TSLM = (mu_A/mu)*(exp(abs(Rx_Compd_TSLM)*log(1+mu)/mu_A)-1).*sign(Rx_Compd_TSLM);
         
         ber_temp_Compd = TSLM_Conv(T_TSLM(:,:,index_min),Rx_DeCompd_TSLM, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
         %ber_temp_Compd = OFDM_Con(SLM_Phase,Rx_DeCompd, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
         ber_count_Compd_TSLM = ber_count_Compd_TSLM+ber_temp_Compd;

    end
    counter=[k_EbN0 length(EbN0_dB)];
    disp(counter);
    % 2 bits per symbol
    ber = ber_count/symbol_total/N_iter/2;
    bers(k_EbN0) = ber;

    ber_Compd_TSLM=ber_count_Compd_TSLM/symbol_total/N_iter/2;
    bers_Compd_TSLM(k_EbN0)= ber_Compd_TSLM;
end

%% 3 Plot
M = N_LEVEL;
Pe_AWGN = ber_QAM(0:EbN0_dB(end),M,'AWGN');
%Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'Rayleigh');
figure(2);
% AWGN Analytic
semilogy(0:EbN0_dB(end),Pe_AWGN,'r-');hold on;grid on;
semilogy(EbN0_dB,bers,'b-*');% Simulation
semilogy(0:EbN0_dB(end),bers_Compd_TSLM, 'p-');% SLM Simulation

legend('AWGN Analytic', 'AWGN Simulation', 'Companding + TSLM Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 15 1e-7 1]);

figure(1);

papr_db=10*log10(Peak./Average);
[f,x]=ecdf(papr_db);
ccdf=1-f;

papr_db_Compd_TSLM=10*log10(Peak_Compd_TSLM./Average_Compd_TSLM);
[f_Compd_TSLM,x_Compd_TSLM]=ecdf(papr_db_Compd_TSLM);
ccdf_Compd_TSLM=1-f_Compd_TSLM;

marker_idx = 1:1e3:length(x);
semilogy(x,ccdf,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Compd_TSLM,ccdf_Compd_TSLM,'*-r','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
%legend('CCDF of Original', 'CCDF of Clipping','CCDF of SLM');
legend('CCDF of Original', 'CCDF of Companding + TSLM');