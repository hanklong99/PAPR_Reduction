%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
k=6;
n_fft=128;
l=k*n_fft*1e3;
EbNodB=0:2:20;
for n=1:length(EbNodB);n
ber(n)=M_QAM_OFDM_fading(l,n_fft,EbNodB(n));
end;
semilogy(EbNodB,ber,'O-');
grid on
xlabel('EbNo')
ylabel('BER')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%