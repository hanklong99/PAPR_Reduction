clear;
close all;
clc

%% Initialization

n_tap=7;

EbN0_dB = 0:1:30;
N_iter = 1e3; % repeat N_iter
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

x_GI = complex(zeros(Nsym * Nframe,1),0);
Y = complex(zeros(symbol_total,1),0);
symbol_output = zeros(symbol_total,1);

ofdm_ray = zeros(Nframe,Nsym);
ofdm_ray_rx= zeros(Nframe,Nsym+n_tap-1);
H=zeros(Nframe, Nfft);

Rx_whole = zeros(1,Nframe*Nsym);

%SLM
N_SLM = 10; % number of sets of SLM Phase
QPSK_SLM=zeros(symbol_total,N_SLM);
bers_SLM_10=zeros(size(EbN0_dB));
Power_SLM=zeros(N_SLM, Nsym*Nframe);
ofdm_SLM=zeros(N_SLM,Nsym*Nframe);
PAPR_SLM=zeros(N_SLM,N_iter);
min_PAPR = zeros(N_iter,1);
PAPR_Total=zeros(N_iter*length(EbN0_dB),1);


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
  
        %Calculate BER
        SLM_Phase=zeros(length(data_QPSK), 1); % Clipping do not need SLM Phase

         % Selective Mapping
         gray_P = randi([0 N_LEVEL-1],length(data_QPSK), N_SLM);
         gray_P=exp(1j * gray_P* pi/2 + 1j * pi/4);
         for i=1:N_SLM
             QPSK_SLM(:,i) = data_QPSK.*gray_P(:,i);
             ofdm_SLM(i,:) = OFDM_Gen(QPSK_SLM(:,i),Ng, Nframe, Nfft, Nsym);
             Power_SLM(i,:) = abs(ofdm_SLM(i,:)).^2;
             PAPR_SLM(i,m)=(max(Power_SLM(i,:))/mean(Power_SLM(i,:)));
         end
         min_PAPR(m)=min(PAPR_SLM(:,m));
         [a,b]=find(PAPR_SLM==min_PAPR(m));
         index_min=a;

         ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(Nframe,n_tap)+1j*randn(Nframe,n_tap));
        
         for i=1:Nframe
         ofdm_ray(i,:)=ofdm_SLM(index_min, (i-1)*Nsym+1:i*Nsym);
         ofdm_ray_rx(i,:) = conv(ht(i,:),ofdm_ray(i,:));
         H(i,:)=fft(ht(i,:),Nfft);
         end
         
        y =  ofdm_ray_rx;

        % AWGN noise
        % sqrt(2) normalization
        complex_noise = sqrt(N0(k_EbN0))*(randn(size(y))+1j*randn(size(y)))/sqrt(2); 
        %Add Noise

         Rx=y+complex_noise;
        
         for i=1:Nframe
             Rx_whole((i-1)*Nsym+1:i*Nsym)=Rx(i,1:Nsym);
         end

        Y = zeros(symbol_total, 1);
        for k = 1:Nframe % N grames
            kk2 = (k-1) * Nsym + (1:Nsym);
            kk1 = (k-1) * Nfft + (1:Nfft);
            Y(kk1) = fft(remove_GI(Ng,Nsym,Rx_whole(kk2)))./H(k,:); 
        end
        
        Y = Y./gray_P(:,index_min);
        
% DeMod
        y_real = real(Y);
        y_imag = imag(Y);

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
% Calculate BER SLM

    ber_SLM=ber_count/symbol_total/N_iter/2;
    bers_SLM_10(k_EbN0)= ber_SLM;
    
    PAPR_Total(1+(k_EbN0-1)*N_iter:k_EbN0*N_iter)=min_PAPR;
    
end

%% 4??Plot
M = N_LEVEL;

Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'Rayleigh');
figure(1);
% Rayleigh Analytic
semilogy(0:EbN0_dB(end),Pe_Ray,'r-');hold on;grid on;

semilogy(0:EbN0_dB(end),bers_SLM_10, 'p-');% SLM Simulation

legend('Rayleigh Analytic', 'SLM Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 30 1e-7 1]);

figure(2);

papr_db_SLM=10*log10(PAPR_Total);
[f_SLM,x_SLM_10]=ecdf(papr_db_SLM);
ccdf_SLM_10=1-f_SLM;

marker_idx = 1:1e3:length(x_SLM_10);
semilogy(x_SLM_10,ccdf_SLM_10,'*-r','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend('CCDF of SLM');