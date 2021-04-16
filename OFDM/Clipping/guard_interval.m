function x_GI = guard_interval( Ng,N_FFT,x )
% 作用：加入CP
% 说明：输入N_FFT个点，输出N_SYM个点
% 版本号：v1.0
% 开始时间：2015年10月12日21:55:21

x_GI = [ x(N_FFT-Ng+1:N_FFT); x];

end