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

bers = zeros(size(EbN0_dB));
bers_clip = zeros(size(EbN0_dB));
Power = zeros(Nsym*Nframe,N_iter);
Average=zeros(N_iter,1);
Peak=zeros(N_iter,1);
Power_clip = zeros(Nsym*Nframe,N_iter);
Average_clip=zeros(N_iter,1);
Peak_clip=zeros(N_iter,1);

%% 1¡¢Calculate the amplitude of Noise
% Alternative method
EbN0=10.^(EbN0_dB/10);
%1/64 is the symbol energy
sigPow = 1/64;
Eb=sigPow/Nbps*T;  %bit energy is half of symbol energy
N0=Eb./EbN0; %sqrt(2) for normalizing the noise magnitude


%% 2¡¢Generate signal again
N_block = N_iter;
%gain = 1/sqrt(2)*[randn(1, N_block) + 1j*randn(1, N_block)];
for k_EbN0 = 1:length(EbN0_dB)
    
    ber_count = 0; % error bit number
    ber_count_clip=0;
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

        % Clipping
        CR=0.001;
        avg=mean(abs(ofdm));
        dev=ofdm-avg;
        sigma=sum(dev*dev')/length(ofdm);
        clipped=ofdm;
        %CL=CR*sigma+avg;
        CL=avg;
        for i=1:length(clipped)
            if abs(clipped(i))>CL
                clipped(i) = (clipped(i)/abs(clipped(i)))*CL;
            end
        end
        Power_clip(:,m)=abs(clipped).^2;
        Average_clip(m)=mean(Power_clip(:,m));
        Peak_clip(m)=max(Power_clip(:,m));
        Rx=clipped+complex_noise;
        ber_temp_clip = OFDM_Con(SLM_Phase, Rx, symbol, Ng, Nframe, Nfft, Nsym, symbol_total);
        ber_count_clip=ber_count_clip+ber_temp_clip;

    end
    counter=[k_EbN0 length(EbN0_dB)];
    disp(counter);
    % 2 bits per symbol
    ber = ber_count/symbol_total/N_iter/2;
    bers(k_EbN0) = ber;

    ber_clip=ber_count_clip/symbol_total/N_iter/2;
    bers_clip(k_EbN0)= ber_clip;
end

%% 4¡¢Plot
M = N_LEVEL;
Pe_AWGN = ber_QAM(0:EbN0_dB(end),M,'AWGN');
%Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'Rayleigh');
figure(2);
% AWGN Analytic
semilogy(0:EbN0_dB(end),Pe_AWGN,'r-');hold on;grid on;
semilogy(EbN0_dB,bers,'b-*');% Simulation
semilogy(0:EbN0_dB(end),bers_clip, 'p-');% Clipping Simulation

%legend('AWGN Analytic', 'AWGN Simulation', 'Clipping Simulation');
legend('AWGN Analytic', 'AWGN Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 15 1e-7 1]);

figure(1);

papr_db=10*log10(Peak./Average);
[f,x]=ecdf(papr_db);
ccdf=1-f;

papr_db_clip=10*log10(Peak_clip./Average_clip);
[f_clip,x_clip]=ecdf(papr_db_clip);
ccdf_clip=1-f_clip;

marker_idx = 1:1e3:length(x);
semilogy(x,ccdf,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_clip,ccdf_clip,'*-g','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
%legend('CCDF of Original', 'CCDF of Clipping','CCDF of SLM');
%legend('CCDF of Original', 'CCDF of Clipping');
legend('CCDF of Original');