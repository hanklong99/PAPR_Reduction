%% Initialization for transmission
    clear all;
    close all;
    clc;
    sigEng = 0;
    N=64;
    K=4;
    Nbps = 2;
    N_iter = 1e4;
    EbN0_dB = 0:1:30;
    Frame=5;
    
    n_tap=7;
    
    Tx_ray = zeros(1,K*N+(Frame*2-1)*N/2+N/2+n_tap-1);
    
    %double sampling
    OQAM_Frame=Frame*2;
    symbol_total = N*Frame;
    % PPN1 and PPN2 and overlapped Frame
    Tx=zeros(1,K*N+(Frame*2-1)*N/2+N/2);

    OQAM=zeros(N,Frame*2); %OQAM signals
                        %sampling by 2*
    data_QAM=zeros(N,Frame);
    symbol=zeros(symbol_total,1);
    bers_fbmc = zeros(1,length(EbN0_dB));
    ber_count=0;
    symbol_output=zeros(symbol_total,1);
    Power = zeros(K*N+(Frame*2-1)*N/2+N/2,N_iter);
    Average=zeros(N_iter,1);
    Peak=zeros(N_iter,1);


%% Prototype Filter (cf M. Bellanger, Phydyas project)
    H1=0.971960;
    H2=sqrt(2)/2;
    H3=0.235147;
    % normalization factor
    factech=1+2*(H1+H2+H3);

    % impulse response
    hef(1:K*N)=0;
    for k=1:K*N-1
    hef(1+k)=1-2*H1*cos(pi*k/(2*N))+2*H2*cos(pi*k/N)-2*H3*cos(pi*k*3/(2*N));
    end

    %normalization
    hef=hef/factech;
    h=hef;  
    figure(1);
    plot(h);
    title ('Phydyas Filter Impulse response');

%% Calculate Aeverage Power
    for m = 1:N_iter
    Tx=zeros(1,K*N+(Frame*2-1)*N/2+N/2);

    % Modulator
        for nframe=1:Frame
    
            symbol(1+(nframe-1)*N:nframe*N) = randi(2^Nbps,N,1) - 1; % information source
            % Info Source & Mod
            graycode = symbol(1+(nframe-1)*N:nframe*N);
            % Gray code (according to website)
            index_3 = find(symbol(1+(nframe-1)*N:nframe*N)==3);
            index_2 = find(symbol(1+(nframe-1)*N:nframe*N)==2);
            graycode(index_3) = 2;
            graycode(index_2) = 3;
            % mapped QAM data
            data_QAM(:,nframe)= exp(1j * graycode* pi/2 + 1j * pi/4);    
    
            % OQAM Modulator
            for k=0:2:N-1 % for k even (k from 0 to M-1)
                OQAM(k+1,2*(nframe-1)+1:2*nframe)=[real(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1))) imag(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1)+1))];
                % frame number is from 0
            end
            % scatterplot(v(:,1));
            for k=1:2:N-1 % for k odd
                OQAM(k+1,2*(nframe-1)+1:2*nframe)=[imag(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1))) real(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1)+1))];
            end
    
        end
    
        for nframe_oqam=1:OQAM_Frame
    
    
            x_fbmc=ifft(OQAM(:,nframe_oqam));
    
            % Duplication of the signal
            x4=[x_fbmc.' x_fbmc.' x_fbmc.' x_fbmc.'];
    
            % Apply filter
            signal=x4.*h;
    
            % Transmitted signal
            % Overlap by N/2
    
            Tx(1+(nframe_oqam-1)*N/2:(nframe_oqam-1)*N/2+K*N)=Tx(1+(nframe_oqam-1)*N/2:(nframe_oqam-1)*N/2+K*N)+signal;
        end
        sigEng = sigEng + sum(abs(Tx).^2);
    end

sigPow = sigEng/N/Frame/N_iter;
EbN0=10.^(EbN0_dB/10);
Eb=sigPow/Nbps;
N0=Eb./EbN0;

%% Full transmission
for k_EbN0 = 1:1:length(EbN0_dB)
    k_EbN0

    ber_count = 0; % error bit number

    for m = 1:N_iter
        Tx=zeros(1,K*N+(Frame*2-1)*N/2+N/2);

%% Modulator
        for nframe=1:Frame

            symbol(1+(nframe-1)*N:nframe*N) = randi(2^Nbps,N,1) - 1; % information source
            % Info Source & Mod
            graycode = symbol(1+(nframe-1)*N:nframe*N);
            % Gray code (according to website)
            index_3 = find(symbol(1+(nframe-1)*N:nframe*N)==3);
            index_2 = find(symbol(1+(nframe-1)*N:nframe*N)==2);
            graycode(index_3) = 2;
            graycode(index_2) = 3;
            % mapped QAM data
            data_QAM(:,nframe)= exp(1j * graycode* pi/2 + 1j * pi/4);    

            % OQAM Modulator
            for k=0:2:N-1 % for k even (k from 0 to M-1)
                OQAM(k+1,2*(nframe-1)+1:2*nframe)=[real(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1))) imag(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1)+1))];
                % frame number is from 0
            end
            % scatterplot(v(:,1));
            for k=1:2:N-1 % for k odd
                OQAM(k+1,2*(nframe-1)+1:2*nframe)=[imag(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1))) real(data_QAM(k+1,nframe))*(j^(k+2*(nframe-1)+1))];
            end

        end

        for nframe_oqam=1:OQAM_Frame


            x_fbmc=ifft(OQAM(:,nframe_oqam));

            % Duplication of the signal
            x4=[x_fbmc.' x_fbmc.' x_fbmc.' x_fbmc.'];

            % Apply filter
            signal=x4.*h;

            % Transmitted signal
            % Overlap by N/2

            Tx(1+(nframe_oqam-1)*N/2:(nframe_oqam-1)*N/2+K*N)=Tx(1+(nframe_oqam-1)*N/2:(nframe_oqam-1)*N/2+K*N)+signal;
        end
        
%% Rayleight Fading Multipath: Filter
        ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(1,n_tap)+1j*randn(1,n_tap));
        H=fft(ht, length(Tx_ray));
        Complex_noise=sqrt(N0(k_EbN0)/2)*(randn(1,length(Tx_ray))+1j*randn(1,length(Tx_ray)));
        Tx_ray=conv(ht,Tx)+Complex_noise;
        Tx_ray=ifft(fft(Tx_ray)./H);

        Power(:,m)=abs(Tx).^2;
        Average(m)=mean(Power(:,m));
        Peak(m)=max(Power(:,m));

%% Demodulator
        out=zeros(N,Frame);
        sest=zeros(N,OQAM_Frame);

        for nframe_oqam=1:OQAM_Frame
            r=Tx_ray(1+(nframe_oqam-1)*N/2:(nframe_oqam-1)*N/2+K*N).*h; % We apply the filter

            u=zeros(1,N);
            for k=1:4
                u=u+r(1,1+(k-1)*N:k*N);                  % Summation
            end

            u=u.';
            sest(:,nframe_oqam)=fft(u);
            sest(:,nframe_oqam)=sest(:,nframe_oqam)/0.6863;                           % Normalization
            
        end
        % from 0: because OQAM-1 will not be able to reach 8
        for nframe_oqam=0:OQAM_Frame-2

            % OQAM demodulation
            if (rem(nframe_oqam,2)==0)
            for k=0:2:N-1 % for k even
                c=[sest(k+1,nframe_oqam+1)*((-1*1i)^(k+nframe_oqam)) sest(k+1,nframe_oqam+2)*((-1*1i)^(k+nframe_oqam+1))];
                % frame number is from 0
                c=real(c);
                out(k+1,nframe_oqam/2+1)=c(1)+1i*c(2);
            end

            for k=1:2:N-1 % for  k odd
                c=[sest(k+1,nframe_oqam+1)*((-1*1i)^(k+nframe_oqam)) sest(k+1,nframe_oqam+2)*((-1*1i)^(k+nframe_oqam+1))];
                c=real(c);
                out(k+1,nframe_oqam/2+1)=c(2)+1i*c(1);
            end
            end

        end  

%% DeMod-QAM
        y_real = real(out);
        y_imag = imag(out);

%% Decide the value
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
%% errors locations
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
    figure(1);
    plot(out,'xk');
    ber = ber_count/symbol_total/N_iter/2;
    bers_fbmc(k_EbN0) = ber;
end



%% Plot



figure(3);
semilogy(0:EbN0_dB(end),bers_fbmc,'p-');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 15 1e-3 1]);
grid on;

figure(4);
papr_db=10*log10(Peak./Average);
[f,x_fbmc]=ecdf(papr_db);
ccdf_fbmc=1-f;
marker_idx = 1:1e3:length(x_fbmc);
semilogy(x_fbmc,ccdf_fbmc,'-s','MarkerIndices',marker_idx);grid on; hold on;
xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend('CCDF of Original');