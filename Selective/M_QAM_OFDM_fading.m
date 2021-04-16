%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO SIMULATE PERFORMANCE OF 64-OFDM IN TIME VARYING FREQUENCY SELECTIVE CHANNEL
% n_bits: Input, length of binary sequence
% n_fft: Input, length of FFT (Fast Fourier Transform)
% EbNodB: Input, energy per bit to noise power spectral density ratio
% ber: Output, bit error rate
% Copyright RAYmaps (www.raymaps.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ber]= M_QAM_OFDM_fading(n_bits,n_fft,EbNodB)

Eb=7;
M=64;
k=log2(M);
n_cyc=32;
EbNo=10^(EbNodB/10);
x=transpose(round(rand(1,n_bits)));
% h1=qammod(M);
% h1.inputtype='bit';
% h1.symbolorder='gray';
y=qammod(x,M);
n_sym=length(y)/n_fft;
n_tap=7;

for n=1:n_sym
    s_ofdm=sqrt(n_fft)*ifft(y((n-1)*n_fft+1:n*n_fft),n_fft);
    %s_ofdm_cyc=[s_ofdm(n_fft-n_cyc+1:n_fft); s_ofdm];
    s_ofdm_cyc=[s_ofdm];
    ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(1,n_tap)+j*randn(1,n_tap));
    Hf=fft(ht,n_fft);
    r_ofdm_cyc=conv(s_ofdm_cyc,ht);
    %r_ofdm_cyc=(r_ofdm_cyc(1:n_fft+n_cyc));
    r_ofdm_cyc=(r_ofdm_cyc(1:n_fft));
    %wn=sqrt((n_fft+n_cyc)/n_fft)*(randn(1,n_fft+n_cyc)+j*randn(1,n_fft+n_cyc));
    wn=sqrt((n_fft+n_cyc)/n_fft)*(randn(1,n_fft)+j*randn(1,n_fft));
    r_ofdm_cyc=r_ofdm_cyc+sqrt(Eb/(2*EbNo))*wn.';
    %r_ofdm=r_ofdm_cyc(n_cyc+1:n_fft+n_cyc);
    r_ofdm=r_ofdm_cyc;
    s_est((n-1)*n_fft+1:n*n_fft)=(fft(r_ofdm,n_fft)/sqrt(n_fft))./Hf.';
end

% h2=qamdemod(M);
% h2.outputtype='bit';
% h2.symbolorder='gray';
% h2.decisiontype='hard decision';
z=qamdemod(s_est.',M);
ber=(n_bits-sum(x==z))/n_bits
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%