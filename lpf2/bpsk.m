

signal = fopen('outputp.txt');
s = textscan(signal,'%f\t%f');
theory = fopen('theoretical_output.txt');
t = textscan(theory,'%f\t%f');

figure();
plot(s{1},s{2}); xlabel('SNR, dB'); ylabel('Error probability');
hold on;
plot(t{1},t{2});
hold off;

fclose(signal);

















































