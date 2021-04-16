function ofdm = OFDM_Gen(data_QPSK, Ng, Nframe, Nfft, Nsym)

% scatterplot(data_QPSK); %constellation map

for k = 1:Nframe % single frame
    % effective symbols
    kk1 = (k-1) * Nfft + (1:Nfft);
    % total symbols
    kk2 = (k-1) * Nsym + (1:Nsym);
    % the data in current frame
    X = data_QPSK(kk1);
    % OFDM signals
    x_ofdm = ifft(X);
    % added guard interval
    ofdm(kk2) = guard_interval(Ng,Nfft,x_ofdm);
end
end