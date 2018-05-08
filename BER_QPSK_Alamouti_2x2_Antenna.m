% Simulation BER of QPSK Alamouti STBC with two transmit and receive antenna
%In this discussion, we will assume that the channel is a flat fading Rayleigh multipath channel and the modulation is QPSK.
%All of This Code was Written by Mostafa Amin-Naji  2017/08/01
% Symorder of simoulation is for 'bin' and 'gray';
% For contact me: Mostafa.Amin.Naji@gmail.com
% My Website: Amin-Naji.com


clc
clear 
close all;
% 
% symorder='gray';  % 'bin' or 'gray'
SNR_Range=0:2.5:12.5;
%bits per symbol
b=2;

for ii=1:length(SNR_Range)
    SNR=SNR_Range(ii);
num_errors=200;

cycle = 0;
error_num_bin=0;

error_num_gray=0;
error_num_normal=0;
h = waitbar(0,'1','Name','Mostafa Amin-Naji  Adv. Comm.');
avrage_p_y11=0;
avrage_p_y21=0;
avrage_p_y12=0;
avrage_p_y22=0;

avrage_p_n1t1=0;
avrage_p_n1t2=0;
avrage_p_n2t1=0;
avrage_p_n2t2=0;
    s=1;
while ((error_num_bin<num_errors) )
s=s+1;



%bits
bitsnumber=4;


%QPSK
M=4;
x = randi([0 M-1],bitsnumber/2,1);  % Random symbols
symbol_bin = pskmod(x,M,pi/4,'bin');
symbol_gray = pskmod(x,M,pi/4,'gray');
symbol_normal=symbol_bin;

%H
H11=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);
H12=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);
H21=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);
H22=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);

%Noise
N11=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);
N12=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);
N21=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);
N22=((randn(1,1))+1i*(randn(1,1)))./sqrt(2);



for z=1:bitsnumber/4
x1_bin=symbol_bin(2*z-1,1);
x2_bin=symbol_bin(2*z,1);
x1_gray=symbol_gray(2*z-1,1);
x2_gray=symbol_gray(2*z,1);

x1_normal=symbol_normal(2*z-1,1);
x2_normal=symbol_normal(2*z,1);


% 

sigma=sqrt(((10^((SNR/10))))/2);
H11(z)=sigma.*H11(z);
H12(z)=sigma.*H12(z);
H21(z)=sigma.*H21(z);
H22(z)=sigma.*H22(z);

Y11_bin=H11.*x1_bin+H12.*x2_bin+N11;
Y12_bin=H21.*x1_bin+H22.*x2_bin+N12;
Y21_bin=H11.*(-conj(x2_bin))+(H12.*conj(x1_bin))+N21;
Y22_bin=H21.*(-conj(x2_bin))+H22.*conj(x1_bin)+N22;

Y11_gray=H11.*x1_gray+H12.*x2_gray+N11;
Y12_gray=H21.*x1_gray+H22.*x2_gray+N12;
Y21_gray=H11.*(-conj(x2_gray))+(H12.*conj(x1_gray))+N21;
Y22_gray=H21.*(-conj(x2_gray))+H22.*conj(x1_gray)+N22;


Y1_normal=H11.*x1_normal+N11;
Y2_normal=H22.*x2_normal+N22;
Y_normal(1,1)=Y1_normal;
Y_normal(2,1)=Y2_normal;
end



%++++++++++++++++++++++++++++++++++++++++


avrage_p_y11=abs(Y11_gray)^2+avrage_p_y11;
avrage_p_y21=abs(Y21_gray)^2+avrage_p_y21;
avrage_p_y12=abs(Y12_gray)^2+avrage_p_y12;
avrage_p_y22=abs(Y22_gray)^2+avrage_p_y22;


avrage_p_n1t1=avrage_p_n1t1+abs(N11)^2;
avrage_p_n1t2=avrage_p_n1t2+abs(N12)^2;
avrage_p_n2t1=avrage_p_n2t1+abs(N21)^2;
avrage_p_n2t2=avrage_p_n2t2+abs(N22)^2;



L=length(Y11_bin)*2;
% 



s1_bin(1,1)=(conj(H11(z)).*Y11_bin(z)+H12(z).*conj(Y21_bin(z))+conj(H21(z)).*(Y12_bin(z))+H22(z).*conj(Y22_bin(z)));
s1_bin(2,1)=(conj(H12(z)).*Y11_bin(1,z)-H11(z).*conj(Y21_bin(z))+conj(H22(z)).*(Y12_bin(z))-H21(z).*conj(Y22_bin(z)));
s1_gray(1,1)=(conj(H11(z)).*Y11_gray(z)+H12(z).*conj(Y21_gray(z))+conj(H21(z)).*(Y12_gray(z))+H22(z).*conj(Y22_gray(z)));
s1_gray(2,1)=(conj(H12(z)).*Y11_gray(1,z)-H11(z).*conj(Y21_gray(z))+conj(H22(z)).*(Y12_gray(z))-H21(z).*conj(Y22_gray(z)));

output_qpsk_alamouti_bin = pskdemod(s1_bin,M,pi/4,'bin');
% output_qpsk_bin = pskdemod(s1_bin,M,pi/4,'bin');


output_qpsk_alamouti_gray = pskdemod(s1_gray,M,pi/4,'gray');
% output_qpsk_gray = pskdemod(s1_gray,M,pi/4,'gray');

output_qpsk_normal = pskdemod(Y_normal,M,pi/4,'bin');

waitbar(error_num_bin / num_errors,h,sprintf('calculating BER for SNR= %d',SNR))
% Check symbol error rate.
[num_bin,~] = biterr(x,output_qpsk_alamouti_bin);
[num_gray,~] = biterr(x,output_qpsk_alamouti_gray);
[num_normal,~] = biterr(x,output_qpsk_normal);

error_num_bin = error_num_bin + num_bin;
error_num_gray = error_num_gray + num_gray;
error_num_normal = error_num_normal + num_normal;

cycle = cycle+4;
end
BER_alamouti_bin=error_num_bin/cycle
BER_alamouti_gray=error_num_gray/cycle
BER_normal=error_num_normal/cycle
delete(h) 
Final_BER_alamouti_bin(ii)=BER_alamouti_bin;
Final_BER_alamouti_gray(ii)=BER_alamouti_gray;
Final_BER_normal(ii)=BER_normal;
avarage_p_y=(avrage_p_y11+avrage_p_y21+avrage_p_y12+avrage_p_y22)/8;
avrage_p_n=(avrage_p_n2t2+avrage_p_n1t1+avrage_p_n2t1+avrage_p_n1t2)/8;

SNR_db(1,ii)=10*log10((avarage_p_y-avrage_p_n)/avrage_p_n)
end
figure (1), clf
semilogy(SNR_Range,Final_BER_alamouti_bin,'b')
hold on
semilogy(SNR_Range,Final_BER_alamouti_gray,'r')
% semilogy(SNR_Range,Final_BER_normal,'g')
axis([SNR_Range([1 end]) 1e-5 1e0]);
grid on;  xlabel('SNR'); ylabel('BER')

% figure
% SNRs=0:2.5:12.5;
%  ber = berfading(SNRs/2,'PSK',4,2);
%     semilogy(SNRs/2,ber);
% bertool
