function [B_P1,B_f] = generate_fft(data,fs)

B_fft = fft(data);
B_P2 = abs(B_fft)/(length(data));
B_P1 = B_P2(1:((length(data))/2) + 1);
B_P1(2:end-1) = 2*B_P1(2:end-1);
B_f = fs*(0:(length(data))/2)/(length(data)-1);

figure;
semilogy(B_f,B_P1);
set(gca,'YScale','log','YMinorTick','on')
    title('Spectra of B Field Magnitude')