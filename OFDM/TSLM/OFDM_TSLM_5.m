clear;
close all;
clc

%% Initialization
EbN0_dB = 0:1:30;
N_iter = 1e3; % repeat N_iter
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
Nframe = 5;       % the number of frames
symbol_total = Nfft * Nframe; % information source

n_tap=7;


%initialize the arrays
x_GI = complex(zeros(Nsym * Nframe,1),0);
Y = complex(zeros(symbol_total,1),0);
symbol_output = zeros(symbol_total,1);

%TSLM
N_TSLM = 5; % number of sets of SLM Phase
bers_TSLM_5=zeros(size(EbN0_dB));
Power_TSLM=zeros(N_TSLM, Nsym*Nframe);
ofdm_TSLM=zeros(N_TSLM,Nsym*Nframe);
PAPR_TSLM=zeros(N_TSLM,N_iter);
min_PAPR = zeros(N_iter,1);
PAPR_Total=zeros(N_iter*length(EbN0_dB),1);
T_TSLM = zeros(Nfft,Nfft, N_TSLM);
F = zeros(Nfft,Nfft, N_TSLM);
s_GI = zeros(Nsym*Nframe,1);

ofdm_ray = zeros(Nframe,Nsym);
ofdm_ray_rx= zeros(Nframe,Nsym+n_tap-1);
H=zeros(Nframe, Nfft);

Rx_whole = zeros(1,Nframe*Nsym);


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
    ber_count_TSLM=0;
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

             for k=1:Nframe
             % effective symbols
             kk1 = (k-1) * Nfft + (1:Nfft);
             % total symbols
             kk2 = (k-1) * Nsym + (1:Nsym);
             t_QPSK = Inv_D*data_QPSK(kk1);
             s = T_TSLM(:,:,i)*t_QPSK;
             s_GI(kk2) = guard_interval(Ng,Nfft,s);
             end
             ofdm_TSLM(i,:)=s_GI;
             Power_TSLM(i,:) = abs(ofdm_TSLM(i,:)).^2;
             PAPR_TSLM(i,m)=(max(Power_TSLM(i,:))/mean(Power_TSLM(i,:)));
         end
         
         min_PAPR(m)=min(PAPR_TSLM(:,m));
         [a,b]=find(PAPR_TSLM==min_PAPR(m));
         index_min=a;


         ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(Nframe,n_tap)+1j*randn(Nframe,n_tap));
        
         for i=1:Nframe
         ofdm_ray(i,:)=ofdm_TSLM(index_min,(i-1)*Nsym+1:i*Nsym);
         ofdm_ray_rx(i,:) = conv(ht(i,:),ofdm_ray(i,:));
         H(i,:)=fft(ht(i,:),Nfft);
         end

        % Channel signal
        y = ofdm_ray_rx;

        %AWGN noise
        % sqrt(2) normalization
        complex_noise = sqrt(N0(k_EbN0))*(randn(size(y))+1j*randn(size(y)))/sqrt(2); 
        %Add Noise


         Rx_TSLM = y + complex_noise;
         
        for i=1:Nframe
            Rx_whole((i-1)*Nsym+1:i*Nsym)=Rx_TSLM(i,1:Nsym);
        end
         
         data_QPSK = zeros(Nfft*Nframe,1);
         TSLM_Phase = T_TSLM(:,:,index_min);
        
         for k = 1:Nframe % N grames
             kk2 = (k-1) * Nsym + (1:Nsym);
             kk1 = (k-1) * Nfft + (1:Nfft);
             
             t_QPSK(kk1) = remove_GI(Ng,Nsym,Rx_whole(kk2));
             t_QPSK(kk1) = TSLM_Phase\t_QPSK(kk1);
             data_QPSK(kk1)= D*t_QPSK(kk1);
             data_QPSK(kk1)=data_QPSK(kk1)./H(k,:).';
         end

% DeMod
         y_real = real(data_QPSK);
         y_imag = imag(data_QPSK);
 
% Decide the value
         for n = 1:symbol_total
             r_value = y_real(n);
             i_value = y_imag(n);
             if r_value > 0 && i_value >0
                 symbol_output(n) = 0;
             elseif r_value < 0 && i_value < 0
                 symbol_output(n) = 3;
             elseif r_value < 0 && i_value > 0
                 symbol_output(n) = 1;
             else
                 symbol_output(n) = 2;
             end
         end
% errors locations
         error_index = find(symbol_output~=symbol);
         
         for error_k = 1:length(error_index)
             k_temp = error_index(error_k);
             error_v = symbol_output(k_temp);
             origin_v = symbol(k_temp);
             
             % 0 --> 1 \ 0 --> 2; 1 bit error
             % 0 --> 3; 2 bits error
             if origin_v == 0
                 if error_v == 1 || error_v == 2
                     ber_count = ber_count + 1;
                 else
                     ber_count = ber_count + 2;
                 end
                 
                 % 1 --> 0 \ 1 --> 3; 1 bit error
                 % 1 --> 2; 2 bits error
             elseif origin_v == 1
                 if error_v == 0 || error_v == 3
                     ber_count = ber_count + 1;
                 else
                     ber_count = ber_count + 2;
                 end
                 
                 
                 % 2 --> 0 \ 2 --> 3; 1 bit error
                 % 2 --> 1; 2 bits error
             elseif origin_v == 2
                 if error_v == 0 || error_v == 3
                     ber_count = ber_count + 1;
                 else
                     ber_count = ber_count + 2;
                 end
                 
                 % 3 --> 1 \ 3 --> 2; 1 bit error
                 % 3 --> 0; 2 bits error
             else
                 if error_v == 1 || error_v == 2
                     ber_count = ber_count + 1;
                 else
                     ber_count = ber_count + 2;
                 end
             end
         end
    end
    counter=[k_EbN0 length(EbN0_dB)];
    disp(counter);
    % 2 bits per symbol

% Calculate BER
    ber_TSLM=ber_count/symbol_total/N_iter/2;
    bers_TSLM_5(k_EbN0)= ber_TSLM;
    
    PAPR_Total(1+(k_EbN0-1)*N_iter:k_EbN0*N_iter)=min_PAPR;
end

%% 3 Plot
M = N_LEVEL;

Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'Rayleigh');
figure(1);
% AWGN Analytic
semilogy(0:EbN0_dB(end),Pe_Ray,'r-');hold on;grid on;
semilogy(0:EbN0_dB(end),bers_TSLM_5, 'p-');% SLM Simulation

legend('Rayleigh Analytic', 'TSLM Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 30 1e-5 1]);

figure(2);

papr_db_TSLM=10*log10(PAPR_Total);
[f_TSLM,x_TSLM_5]=ecdf(papr_db_TSLM);
ccdf_TSLM_5=1-f_TSLM;

marker_idx = 1:1e3:length(x_TSLM_5);
semilogy(x_TSLM_5,ccdf_TSLM_5,'*-r','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend( 'CCDF of TSLM');