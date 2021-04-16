function ber_count= OFDM_Con(SLM_Phase, Rx, symbol, Ng, Nframe, Nfft, Nsym, symbol_total)
        symbol_output = zeros(symbol_total,1);  
        % Rx signal
        ber_count=0;
        Y = zeros(symbol_total, 1);
        for k = 1:Nframe % N grames
            kk2 = (k-1) * Nsym + (1:Nsym);
            kk1 = (k-1) * Nfft + (1:Nfft);
            
            Y(kk1) = fft(remove_GI(Ng,Nsym,Rx(kk2)));
        end
        
        if(SLM_Phase ~= 0)
        Y = Y./SLM_Phase;
        end
        % DeMod
        y_real = real(Y);
        y_imag = imag(Y);

        % Decide the value
        for n = 1:symbol_total
            r_value = y_real(n);
            i_value = y_imag(n);
            if r_value > 0 && i_value >0
                symbol_output(n) = 0;
            elseif r_value < 0 && i_value < 0
                symbol_output(n) = 3;
            elseif r_value < 0 && i_value > 0
                symbol_output(n) = 1;
            else
                symbol_output(n) = 2;
            end
        end
        % errors locations
        error_index = find(symbol_output~=symbol);
        
        for error_k = 1:length(error_index)
            k_temp = error_index(error_k);
            error_v = symbol_output(k_temp);
            origin_v = symbol(k_temp);
            
            % 0 --> 1 \ 0 --> 2; 1 bit error
            % 0 --> 3; 2 bits error
            if origin_v == 0
                if error_v == 1 || error_v == 2
                    ber_count = ber_count + 1;
                else
                    ber_count = ber_count + 2;
                end
                
                % 1 --> 0 \ 1 --> 3; 1 bit error
                % 1 --> 2; 2 bits error
            elseif origin_v == 1
                if error_v == 0 || error_v == 3
                    ber_count = ber_count + 1;
                else
                    ber_count = ber_count + 2;
                end
                
                
                % 2 --> 0 \ 2 --> 3; 1 bit error
                % 2 --> 1; 2 bits error
            elseif origin_v == 2
                if error_v == 0 || error_v == 3
                    ber_count = ber_count + 1;
                else
                    ber_count = ber_count + 2;
                end
                
                % 3 --> 1 \ 3 --> 2; 1 bit error
                % 3 --> 0; 2 bits error
            else
                if error_v == 1 || error_v == 2
                    ber_count = ber_count + 1;
                else
                    ber_count = ber_count + 2;
                end
            end
            
        end
        
end