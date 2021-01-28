clc;
close all;
clear all;
 
% Message signal 
L1=44000;
bits=11*L1;
%Generate random message bits
X=randi([0,1],1,bits);  
 
%Generator matrix for hamming code
N=15;
K=11;
%Parity Matrix  
P=[1 1 1 1; 0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0; 0 0 1 1; 0 1 0 1; 0 1 1 0; 1 0 1 0; 1 0 0 1; 1 1 0 0];   
Ik=eye(K);  
Im=eye(N-K);
G=[P Ik];    %Generator matrix
H=[Im P'];   %Parity check matrix

%Error pattern
E=vertcat(zeros(1,15),eye(15));  
S=mod(E*H',2);   %Syndrome
 
% Encoding using Hamming code
 
sig=reshape(X,11,L1);  %11 per column
%Encode column by column
en_sig1=mod(sig'*G,2)'; 
% Encoded signal
en_sig=2*reshape(en_sig1,1,15*L1)-1;  
 
 
%QPSK modulation
Y=en_sig;  %Encoded QPSK input signal
REAL=Y(1:2:end);
IMAG=Y(2:2:end);
U=REAL+1i*IMAG;  %QPSK Modulation

%MRC diversity technique
H1=sqrt(1/2)*(randn(1,1000)+1i*randn(1,1000)); 
H2=sqrt(1/2)*(randn(1,1000)+1i*randn(1,1000));

FINALRECVECTOR=zeros(1,660000);

%QPSK transmitted signal

for i=0:999
  Y1(1,i*330+1:(i+1)*330)=H1(i+1)*U(1,i*330+1:(i+1)*330);   
Y2(1,i*330+1:(i+1)*330)=H2(i+1)*U(1,i*330+1:(i+1)*330);
end

%Generate AWGN for transmitted signal
AW_noise1=randn(1,size(U,2))+1i*randn(1,size(U,2));  
AW_noise2=randn(1,size(U,2))+1i*randn(1,size(U,2));  

for k=1:15
    SNR=10^(k/10);  %dB to linear
    N0=1/SNR;
        
%Add AWGN
COPY1=Y1+sqrt((N0*(15/11))/2)*AW_noise1;  
COPY2=Y2+sqrt((N0*(15/11))/2)*AW_noise2;
        
for i=0:999
  F1(1,i*330+1:(i+1)*330)=conj(H1(i+1))*COPY1(1,i*330+1:(i+1)*330);    F2(1,i*330+1:(i+1)*330)=conj(H2(i+1))*COPY2(1,i*330+1:(i+1)*330);
end

FINALREC=F1+F2;
REALFINAL=real(FINALREC);
IMAGFINAL=imag(FINALREC);
 
%QPSK demodulation
 
for l=1:330000
 
if REALFINAL(1,l)>0
REALFINAL(1,l)=1;
else 
    REALFINAL(1,l)=0;
end
 
if  IMAGFINAL(1,l)>0
IMAGFINAL(1,l)=1;
else 
    IMAGFINAL(1,l)=0;
end
 
end
FINALRECVECTOR(1:2:end)=REALFINAL;
FINALRECVECTOR(2:2:end)=IMAGFINAL;
 
%Decoding using Syndrome
 
R=reshape(FINALRECVECTOR,15,L1).';

%Generate error syndromes
er_syn=mod(R*H',2);  
for k1=1:L1
    for k2=1:size(E,1)
        if S(k2,:)== er_syn(k1,:)
     idxe=k2; %Find the syndrome index 
        end
    end

%Look up the error pattern
    error=E(idxe,:);  
%Error correction
    corr_sig(k1,:)=xor(R(k1,:),error);  
%Keep the message bits
   sig_bit=corr_sig(:,5:15);  
end
 
msg=reshape(sig_bit',1,bits);   
BER(k)=sum(xor(msg,X));
%Calculation of Simulated Bit Error Rate
estimated_BER(k)=BER(k)/bits;   
 
end
k=1:15;

% BER Vs Eb/N0 curve
semilogy(k,estimated_BER,'bs-','LineWidth',2);   
grid on;
xlabel('Eb/N0(dB)');
ylabel('BER');
title('BER Vs Eb/N0, Hamming Code(15,11) coded QPSK');


