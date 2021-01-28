clc;
clear all;
close all; 

%Generate message bits
X=randi([0,1],1,484000);  
Y=(2*X)-1;      %NRZ signal convertion

REAL=Y(1:2:end);
IMAG=Y(2:2:end);
U=REAL+1i*IMAG;      %QPSK Modulation

%MRC diversity technique
H1=sqrt(1/2)*(randn(1,1000)+1i*randn(1,1000));      
H2=sqrt(1/2)*(randn(1,1000)+1i*randn(1,1000));

FINALRECVECTOR=zeros(1,484000);

for i=0:999
%QPSK transmitted signal
Y1(1,i*242+1:(i+1)*242)=H1(i+1)*U(1,i*242+1:(i+1)*242);    
Y2(1,i*242+1:(i+1)*242)=H2(i+1)*U(1,i*242+1:(i+1)*242);
end

for k=1:15
    No(k)=10^(-k/10);  %dB to linear

%Generate AWGN for transmitted signal     
N1=sqrt(No(k)/2)*(randn(1,size(U,2))+1i*randn(1,size(U,2)));        
N2=sqrt(No(k)/2)*(randn(1,size(U,2))+1i*randn(1,size(U,2)));

%Add AWGN
        COPY1=Y1+N1; 
        COPY2=Y2+N2;
 
for i=0:999  F1(1,i*242+1:(i+1)*242)=conj(H1(i+1))*COPY1(1,i*242+1:(i+1)*242);  
F2(1,i*242+1:(i+1)*242)=conj(H2(i+1))*COPY2(1,i*242+1:(i+1)*242);
end

FINALREC=F1+F2;
REALFINAL=real(FINALREC);
IMAGFINAL=imag(FINALREC);
 
%QPSK Demodulation
 
for l=1:242000
 
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
M=mod(FINALRECVECTOR+X,2);
ERRORS=sum(M);
%Calculation of Simulated Bit Error Rate
Pe(k)=ERRORS/484000;  
end

k=1:15;
% BER Vs Eb/N0 curve
semilogy(k,Pe,'rs-','LineWidth',2)  
grid on;
legend('Coded', 'Uncoded'); 
xlabel('Eb/N0(dB)');
ylabel('BER');
title('BER Vs Eb/N0: uncoded QPSK');
