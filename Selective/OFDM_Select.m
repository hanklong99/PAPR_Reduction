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
Ng = (1/4)*Nfft;      % Guard interval
Nsym = Nfft + Ng; % total number of symbols per frame
Nframe = 5;       % the number of frames
symbol_total = Nfft * Nframe; % information source
% Selective channel
% PDP= [0.6 0.3 0.1]; %channel with three taps
% L=Nframe; %number of channel realizations, change this to 1 if only one realization is needed
% N = length(PDP); %number of channel taps
% a = sqrt(PDP); %path gains (tap coefficients a[n])

n_tap=7;


%initialize the arrays
x_GI = complex(zeros(Nsym * Nframe,1),0);
Y = complex(zeros(symbol_total,1),0);
symbol_output = zeros(symbol_total,1);
bers = zeros(size(EbN0_dB));
Power = zeros(Nsym*Nframe,N_iter);
Average=zeros(N_iter,1);
Peak=zeros(N_iter,1);
Power_clip = zeros(Nsym*Nframe,N_iter);
Average_clip=zeros(N_iter,1);
Peak_clip=zeros(N_iter,1);

ofdm_ray = zeros(Nframe,Nsym);
ofdm_ray_rx= zeros(Nframe,Nsym+n_tap-1);
H=zeros(Nframe, Nfft);

Rx_whole = zeros(1,Nframe*Nsym);


%% 1 Calculate the signal power
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
    ofdm=OFDM_Gen(data_QPSK, Ng, Nframe, Nfft, Nsym);
    
    Power(:,m)=abs(ofdm).^2;
    % total energy
    sigEng = sigEng + sum(abs(ofdm).^2); % conjugate to calculate energy
    %Power(:,m)=abs(x_GI).^2;
    Average(m)=mean(Power(:,m));
    Peak(m)=max(Power(:,m));
end

% signal power of one symbol
sigPow = sigEng/Nsym/Nframe/N_iter;

%% 2 Calculate the amplitude of Noise
% Alternative method
EbN0=10.^(EbN0_dB/10);
%1/64 is the symbol energy
Eb=sigPow/Nbps;        %bit energy is half of symbol energy
N0=Eb./EbN0; %sqrt(2) for normalizing the noise magnitude



%% 3 Generate signal again
N_block = N_iter;

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
        ofdm=OFDM_Gen(data_QPSK,Ng, Nframe, Nfft, Nsym);
% Frequency Selective        
        
        %Rayleigh a set of random variables R0,R1,R2,…R{N-1} with unit average power
        %for each channel realization (block type fading)
%         R =1/sqrt(2).*((randn(L,N))+1i*randn(L,N));% + i*(ones(L,N))); % normalized the Rayleigh fading variables
%         taps= repmat(a,L,1).*R; %combined tap weights = a[n]*R[n]
%         h= 1/sqrt(sum(PDP))*taps;%Normalized taps for output with unit average power
        
        ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(Nframe,n_tap)+1j*randn(Nframe,n_tap));
        
        for i=1:Nframe
        ofdm_ray(i,:)=ofdm((i-1)*Nsym+1:i*Nsym);
        ofdm_ray_rx(i,:) = conv(ht(i,:),ofdm_ray(i,:));
        H(i,:)=fft(ht(i,:),Nfft);
        end
        
        %Now, you can use the rest of the code as
        % Channel signal
        y = ofdm_ray_rx;
        %AWGN noise
        complex_noise = sqrt(N0(k_EbN0))*(randn(size(y))+1j*randn(size(y)))/sqrt(2); % sqrt(2) normalization
        
        Rx=y+complex_noise;
        
        for i=1:Nframe
            Rx_whole((i-1)*Nsym+1:i*Nsym)=Rx(i,1:Nsym);
        end
        %SLM_Phase = 0;
        %Y = OFDM_Con(SLM_Phase,Rx_whole, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
        Y = zeros(symbol_total, 1);
        for k = 1:Nframe % N grames
            kk2 = (k-1) * Nsym + (1:Nsym);
            kk1 = (k-1) * Nfft + (1:Nfft);
            Y(kk1) = fft(remove_GI(Ng,Nsym,Rx_whole(kk2)))./H(k,:); 
        end
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
    ber = ber_count/symbol_total/N_iter/2;
    bers(k_EbN0) = ber;
end

%% 4 Plot
M = N_LEVEL;
Pe_AWGN = ber_QAM(0:EbN0_dB(end),M,'Rayleigh');
Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'AWGN');
figure(2);
semilogy(0:EbN0_dB(end),Pe_AWGN,'r-');hold on;grid on;
semilogy(0:EbN0_dB(end),Pe_Ray, 'g-');
semilogy(EbN0_dB,bers,'b-*'),
legend('Rayleigh Analytic', 'AWGN Analytic', 'Ray Select Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 30 1e-7 1]);

% figure(1);
% 
% papr_db=10*log10(Peak./Average);
% [f,x]=ecdf(papr_db);
% ccdf=1-f;
% 
% marker_idx = 1:1e3:length(x);
% semilogy(x,ccdf,'-s','MarkerIndices',marker_idx);grid on; hold on;
% 
% xlabel('PAPR/dB'); ylabel('CCDF');
% xlim([0 13]); ylim([10e-4 1]);
% legend('CCDF of Original');