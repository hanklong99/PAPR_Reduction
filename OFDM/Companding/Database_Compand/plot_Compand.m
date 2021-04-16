figure(1);

semilogy(0:EbN0_dB(end),bers, 'r-<'); hold on; grid on;
semilogy(0:EbN0_dB(end),bers_Compd_2, 'y->'); hold on; grid on;
semilogy(0:EbN0_dB(end),bers_Compd_3, 'm-<');
semilogy(0:EbN0_dB(end),bers_Compd_5, 'k--x');
legend('No Companding','Companding 2', 'Companding 3', 'Companding 5');
xlabel('Eb/N0 [dB]'), ylabel('BER'); axis([0 30 1e-5 1]);


figure(2);
marker_idx = 1:1e3:length(x_Compd_2);

semilogy(x,ccdf,'*-r','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Compd_2,ccdf_Compd_2,'*-g','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Compd_3,ccdf_Compd_3,'>-r','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Compd_5,ccdf_Compd_5,'-s','MarkerIndices',marker_idx);grid on; hold on;
xlabel('PAPR/dB'); ylabel('CCDF');
axis([4 11 1e-3 1]);
% xlim([2 11]); ylim([10e-4 1]);
legend('No Companding','Companding 2', 'Companding 3', 'Companding 5');