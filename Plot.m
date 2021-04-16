figure(1);

hold on;
grid on;
% semilogy(0:EbN0_dB(end),Pe_AWGN,'r-+');hold on;grid on;
% semilogy(0:EbN0_dB(end),Pe_Ray, 'g-s');
semilogy(EbN0_dB,bers,'b-*'),
semilogy(0:EbN0_dB(end),bers_clip, 'k--x');
% semilogy(0:EbN0_dB(end),bers_SLM, 'y->');
% semilogy(0:EbN0_dB(end),bers_TSLM, 'm-<');
%semilogy(0:EbN0_dB(end),bers_Compd_TSLM, 'k-p');
semilogy(0:EbN0_dB(end),bers_Compd, 'c-h');
xlim([0 15]); ylim([10e-3 1]);
%legend('AWGN Analytic', 'Rayleigh Analytic', 'AWGN Simulation', 'Clipping Simulation', 'SLM Simulation','TSLM Simulation','Companding&TSLM Simulation','Companding Simulation');
legend( 'Multipath Rayleigh Simulation', 'Clipping Simulation', 'Companding Simulation');
xlabel('E_b/N_0 [dB]'), ylabel('BER'); 


figure(2);
marker_idx = 1:1e3:length(x);
semilogy(x,ccdf,'-s','MarkerIndices',marker_idx);grid on; hold on;
semilogy(x_clip,ccdf_clip,'*-g','MarkerIndices',marker_idx);
% semilogy(x_SLM,ccdf_SLM,'>-r','MarkerIndices',marker_idx);
% semilogy(x_TSLM,ccdf_TSLM,'<-m','MarkerIndices',marker_idx);
%semilogy(x_Compd_TSLM,ccdf_Compd_TSLM,'p-c','MarkerIndices',marker_idx);
semilogy(x_Compd,ccdf_Compd,'h-y','MarkerIndices',marker_idx);

xlabel('PAPR/dB'); ylabel('CCDF');
xlim([0 13]); ylim([10e-4 1]);
legend('CCDF of Original', 'CCDF of Clipping', 'CCDF of Companding');
%legend('CCDF of Original', 'CCDF of Clipping','CCDF of SLM','CCDF of TSLM','CCDF of Companding & TSLM', 'CCDF of Companding');