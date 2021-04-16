figure(1);

semilogy(0:EbN0_dB(end),bers, 'k--x'); hold on; grid on;
semilogy(0:EbN0_dB(end),bers_fbmc, 'r--x'); hold on; grid on;
semilogy(0:EbN0_dB(end),bers_Hybrid_ofdm, 'r->');
semilogy(0:EbN0_dB(end),bers_Hybrid_fbmc, 'g->');
legend('No Reduction OFDM', 'No Reduction FBMC', 'OFDM Hybrid', 'FBMC Hybrid');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); axis([0 30 1e-5 1]);

figure(2);
marker_idx = 1:1e3:length(x);
semilogy(x,ccdf,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_fbmc,ccdf_fbmc,'-r','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Hybrid_ofdm,ccdf_Hybrid_ofdm,'*-g','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_Hybrid_fbmc,ccdf_Hybrid_fbmc,'*-r','MarkerIndices',marker_idx);grid on; hold on;

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend('No Reduction OFDM', 'No Reduction FBMC', 'Hybrid Reduction OFDM', 'Hybrid Reduction FBMC');