clear;
close all;
clc

%% Initialization
n_tap=7;

EbN0_dB = 0:1:100;
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

x_GI = complex(zeros(Nsym * Nframe,1),0);
Y = complex(zeros(symbol_total,1),0);
symbol_output = zeros(symbol_total,1);

ofdm_ray = zeros(Nframe,Nsym);
ofdm_ray_rx= zeros(Nframe,Nsym+n_tap-1);
H=zeros(Nframe, Nfft);

Rx_whole = zeros(1,Nframe*Nsym);


bers_clip = zeros(size(EbN0_dB));
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


        % Clipping
        CR=0.00001;
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
        
        % Pass Frequency Selective Channel
        ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(Nframe,n_tap)+1j*randn(Nframe,n_tap));
        
        for i=1:Nframe
        ofdm_ray(i,:)=clipped((i-1)*Nsym+1:i*Nsym);
        %ofdm_ray(i,:)=ofdm((i-1)*Nsym+1:i*Nsym);
        ofdm_ray_rx(i,:) = conv(ht(i,:),ofdm_ray(i,:));
        H(i,:)=fft(ht(i,:),Nfft);
        end

        % Channel signal
        y = ofdm_ray_rx;
        
        % sqrt(2) normalization
        complex_noise = sqrt(N0(k_EbN0))*(randn(size(y))+1j*randn(size(y)))/sqrt(2); 
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

    ber_clip=ber_count/symbol_total/N_iter/2;
    bers_clip(k_EbN0)= ber_clip;
end

%% 4¡¢Plot

%Pe_Ray = ber_QAM(0:EbN0_dB(end), M, 'Rayleigh');
figure(1);
semilogy(0:EbN0_dB(end),bers_clip, 'p-');% Clipping Simulation

legend('Clip Simulation');
xlabel('Eb/N0 [dB]'), ylabel('BER'); axis([0 30 1e-5 1]);

figure(2);

papr_db_clip=10*log10(Peak_clip./Average_clip);
[f_clip,x_clip]=ecdf(papr_db_clip);
ccdf_clip=1-f_clip;

marker_idx = 1:1e3:length(x_clip);

semilogy(x_clip,ccdf_clip,'*-g','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend('CCDF of Clipping');