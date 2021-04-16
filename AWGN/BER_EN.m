EbNo_dB = 0:0.01:15; %Set the Energy/Power Ratio
EbNo = 10.^(EbNo_dB/10);
BER_U_NRZ = Q(sqrt(EbNo));
BER_B_NRZ = Q(sqrt(2*EbNo));
semilogy(EbNo_dB, BER_U_NRZ);
hold on;
grid on;
semilogy(EbNo_dB, BER_B_NRZ);
legend('Bipolor NRZ','Unipolar NRZ');
xlabel('Eb/No(dB)');
ylabel('BER');