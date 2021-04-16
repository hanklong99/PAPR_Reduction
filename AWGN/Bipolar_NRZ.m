%UnipolarNRZ_AWGN.m Simulate the performance of Unipolar NRZ in AWGN channels
clear all;
close all;
%Eb/N0 in dB
Eb_N0_dB=0:1:10;
%Convert Eb/N0 into the linear scale
Eb_N0=10.^(Eb_N0_dB/10);
%Fix the average bit engergy to be 1
Eb=1;
%Fix the bit duration to be 1
T=1;
%Calculate the signal amplitude
%Or you can also change the amplitude from 2Eb to Eb
A=sqrt(Eb/T);
%Engergy of bit 0
%E0=0;
E0=A^2*T;
%Engergy of bit 1
E1=A^2*T;
%Optimal threshold
Thr=(E1-E0)/2;
%Matched filter impulse response--constant with amplitude A over one bit
%duration.
h=A;
%Calculate the correspnding noise power spectral density
%Bipolar Eb = 2 Unipolar Eb
N0=Eb./Eb_N0;
%Number of different SNRs
len_EbN0=length(Eb_N0);
%Total number of bits per block
N_bit=100000;
%Maximum number of blocks to simulate
N_block=1000;
%Eb/N0 index pointer
EbN0_pointer=1;
temp_EbN0_pointer=EbN0_pointer;
%Number of errors counted for each Eb/N0
errs=zeros(1,len_EbN0);
%Number of blocks simulated for each Eb/N0
block_count=zeros(1,len_EbN0);
%While the Eb/N0 index pointer has not reached the last value, and the
%number of blocks has not exceeded the maximum number, N_block,
%do the iterations
while (EbN0_pointer <= len_EbN0) && (block_count(len_EbN0) < N_block)
 %Generate a binary bit sequence
 D=round(rand(1,N_bit));
 %D=rand(1,N_bit);
 F = zeros(1, N_bit);
 F (D == 1) = 1;
 F (D == 0) = -1;
 
 %Transmitted signal after unipolar NRZ line coding.
 %1 is mapped to amptitude A; 0 is mapped to amptitude 0.
 Tx_data = A*F;
%  Tx_data = A*randsrc(1, N_bit);

 %AWGN with normalised power 1
 Noise=randn(1, N_bit);

 %Simulate different Eb/N0 values
 for n = EbN0_pointer : len_EbN0
 %Standard deviation of AWGN
 sigma_AWGN=sqrt(N0(n)/2);

 %Received signal. The noise power is N0(n).
 Rx_data = Tx_data + sigma_AWGN*Noise;

 %Matched filter output signal
 V=h*Rx_data;

 %Recover the transmit data from the received data.
 %When the MF output is V>Th, output 1; otherwise output 0.
 Recov_data = zeros (1,N_bit);
 Recov_data (V>Thr)=1;
%  Recov_data = sign(Rx_data);


 %Count the number of errors
 errs(n)= errs(n)+sum(abs(Recov_data-D));
%  errs(n) = errs(n) + sum(abs(Recov_data - Tx_data));
 %If more than 500 errors have been counted, move the Eb/N0
 %index pointer to the next Eb/N0
 if errs(n)>=500
 temp_EbN0_pointer = temp_EbN0_pointer+1;
 end
 %Update the nubmer of blocks simulated
 block_count(n)=block_count(n)+1;
 end

 %Update Eb/N0 pointer
 EbN0_pointer=temp_EbN0_pointer;

 block_count
end

%Calculate the numerical BERs for different Eb/N0's. Each block has N_bit bits.
Num_BER = errs./(N_bit*block_count);
%Calculate the analytical BERs for different Eb/N0's
Ana_BER=Q(sqrt(2*Eb_N0));
figure;
semilogy(Eb_N0_dB, Num_BER, '-s');
hold on;
semilogy(Eb_N0_dB, Ana_BER, 'r-*');
grid on;
legend('Bipolar Numerical BER', 'Bipolar Analytical BER');
title('NRZ in AWGN Channels');
xlabel('Eb/N0 (dB)');
ylabel('BER');