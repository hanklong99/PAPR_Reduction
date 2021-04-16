figure(1);

semilogy(0:EbN0_dB(end),bers_TSLM_1, 'k--x'); hold on; grid on;
semilogy(0:EbN0_dB(end),bers_TSLM_5, 'y->');
semilogy(0:EbN0_dB(end),bers_TSLM_10, 'm-<');
semilogy(0:EbN0_dB(end),bers_TSLM_20, 'k-p');
semilogy(0:EbN0_dB(end),bers_TSLM_30, 'c-h');
legend('TSLM_1', 'TSLM_5', 'TSLM_10', 'TSLM_20', 'TSLM_30');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 30 1e-5 1]);

figure(2);
marker_idx = 1:1e3:length(x_TSLM_1);
semilogy(x_TSLM_1,ccdf_TSLM_1,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_TSLM_5,ccdf_TSLM_5,'*-g','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_TSLM_10,ccdf_TSLM_10,'>-r','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_TSLM_20,ccdf_TSLM_20,'<-m','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_TSLM_30,ccdf_TSLM_30,'h-y','MarkerIndices',marker_idx);grid on; hold on;

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend('TSLM_1', 'TSLM_5','TSLM_10','TSLM_20', 'TSLM_30');