figure(1);

semilogy(0:EbN0_dB(end),bers_Compd, 'k--x'); hold on; grid on;
semilogy(0:EbN0_dB(end),bers_Compd_2, 'y->');
semilogy(0:EbN0_dB(end),bers_Compd_3, 'm-<');
legend('Companding 5', 'Companding 2', 'Companding 3');
xlabel('Eb/N0 [dB]'), ylabel('BER'); axis([0 30 1e-5 1]);

figure(2);
marker_idx = 1:1e3:length(x_Compd);
semilogy(x_Compd,ccdf_Compd,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Compd_2,ccdf_Compd_2,'*-g','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Compd_3,ccdf_Compd_3,'>-r','MarkerIndices',marker_idx);grid on; hold on;

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([10 13]); ylim([10e-4 1]);
legend('Companding 5', 'Companding 2', 'Companding 3');