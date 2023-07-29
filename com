Exp--------111111111111111111111111

A = input("Enter the amplitude : ");
F = input("Enter the frequency : ");
t = 0:0.1:10;
%Sine Wave
y1 = A*sin(2*pi*F*t);
subplot(2,2,1);
plot(t,y1);
title("Sine wave");
xlabel("Time");
ylabel("Amplitude");

%Cos wave
y2 = A*cos(2*pi*F*t);
subplot(2,2,2);
plot(t,y2);
title("Cos wave");
xlabel("Time");
ylabel("Amplitude");
%Unit impulse
n =-10:10;
u = [zeros(1,10) 1 zeros(1,10)];
subplot(2,2,1);
plot(n,u);
title("Unit Impulse Wave form in continuse signal");
xlabel("n");
ylabel("u(n)");
%discrete signal
subplot(2,2,3);
plot(n,u);
stem(n,u);
title("Unit Impulse Wave form in discrete signal");
xlabel("n");
ylabel("u(n)");

%Unit step wave
n = -10:10;
u = [zeros(1,10) ones(1,11)];
stem(n,u);
axis([-12 12 -4 4]);
% Square wave 

t=0:0.01:1;
a=2;
b=a*square(2*pi*2*t);
subplot(2,2,1);
plot(t,b);
xlabel('time');
ylabel('Amplitude');
title ('square wave in continous form');
subplot(2,2,2);
stem(t,b);
xlabel('time');
ylabel('Amplitude');
title ('square wave in Discrete form');
 
%Exponential waveform
A = input("Enter the amplitude : ");
F = input("Enter the frequency : ");
t = 0:0.1:10;
%exp Wave
y1 = A*exp(2*pi*F*t);
subplot(2,2,1);
plot(t,y1);
title("EXP wave");
xlabel("Time");
ylabel("Amplitude");




exp2222222222222222222222222222222222
clc;
clear all;
close all;
%using predefined input => x=[0 1 2 3]; 
%x=[0 1 2 3];
% This is for input from the useer => x =input("Enter the sequence : ");
x =input("Enter the sequence : ");
N=length(x);
X=zeros(N,1);
for k=0:N-1
    for n=0:N-1
        X(k+1)=X(k+1)+x(n+1)*exp(-j*2*pi*k*n/N);
    end
end
n=0:N-1;
subplot(2,2,1);
stem(n,x);
title("Input Discrete signal");
xlabel("Time");
ylabel("Amplitude");
% it is the Discrete signal 
grid on;
%Magnitide and phase is the output of DFT process
subplot(2,2,2);
stem(n,abs(X));
title("Magnitude Response for DFT signal");
xlabel("Magnitude");
ylabel("Amplitude");
grid on;
subplot(2,2,3);
stem(n,angle(X));
title("Phase Response for DFT Signal");
xlabel("Time");
ylabel("Phase");
grid on;
y=zeros(N,1);
for n=0:N-1
    for k=0:N-1
        y(n+1)=y(n+1)+X(k+1)*1/N*exp(j*2*pi*k*n/N);
    end
end
disp("IDFT is: ");
disp(y);
k=0:N-1;
subplot(2,2,4);
stem(k,abs(y));
title("IDFT");
xlabel("Time");
ylabel("Amplitude");
%IDFT signal



Ep3 convulationnnnnnnnnnnnnnnnnnnnnn

clc;
x1 = input('Enter the first sequence: ');
subplot(3,1,1);
stem(x1);
ylabel('Amplitube');
title('Plot the first sequence');
x2 = input('Enter the Second sequence: ');
subplot(3,1,2);
stem(x2);
ylabel('Amplitube');
title('Plot the second sequence');
f = conv(x1,x2);
disp(f);
xlabel('time index n');
ylabel('amplitude f');
subplot(3,1,3);
stem(f);
title('linear conv of sequence');








Exp444444444 autocorelation
clc;
clear;
x = input('Enter first ');
subplot(2,1,1);
stem(x);
ylabel('Amplitude');
title('Sequence before auto corelation');
subplot(2,1,2);
y = xcorr(x,x);
stem(y);
ylabel('Amplitude');
title('Sequence after auto corelation');
disp(y);










Exp555555555 fft fft 
clc;
clear all;
close all;
x =input("Enter the sequence: ");
N=length(x);
F=fft(x,N);
n=0:N-1
subplot(2,2,1);
stem(n,x);
xlabel("Time");
ylabel("Amplitide");
title("Input Sequence");
subplot(2,2,2);
stem(n,abs(F));
xlabel("Time");
ylabel("Magnitude");
title(" Magnitude figure");
subplot(2,2,3);
stem(n,angle(F));
xlabel("Time");
ylabel("Phase");
title(" Phase figure");



exp6666666666666 inerpolation

clc;
clear all;
close all;
F=input("Enter the frequency of the signal: ");
P=input("Enter the interpulator factor: ");
N=input("Enter the length of the signal: ");
t=0:1:N-1
%X = sin(2*pi*F*t);
X = sin(2*3.14*F*t);
i=interp(X,P);
subplot(2,1,1);
stem(X);
title("Original signal");
subplot(2,1,2);
stem(i);
title("Interpolated signal");

exp777777777  NNn dft DFT

close all;
clear all;
clc;
x1 = input('Enter the sequence: ');
n = input('Enter the length ');
m = abs(fft(x1,n));
%m =(fft(x1,n));
disp('N Point DFT of a given sequence is : ');
disp(m);
N = 0:1:n-1;
%magnitude plot
subplot(2,2,1);
stem(N,m);
xlabel('Length');
ylabel('Amplitude of x(k)');
title('Magnitude spectrun');
hold on;
%phase plot
an = angle(fft(x1,n));
disp(an);
subplot(2,2,2);
stem(N,an);
xlabel('Length');
ylabel('Phase of x(k)');
title('phase spectrun');



exp888 samplin sampling

clc;
close all;
clear all;
t=-10:.01:10;
T=4;
fm=1/T;
x=cos(2*pi*fm*t);
subplot(2,2,1);
plot(t,x);
xlabel("time");
ylabel("Amplitude");
grid on;
title("Input signal");
n1=-4:1:4;
fs1=1.6*fm;
fs2=2*fm;
fs3=8*fm;
x1=cos(2*pi*fm/fs1*n1);
subplot(2,2,2);
stem(n1,x1);
xlabel("Number of samples");
ylabel("Amplitude");
hold on;
subplot(2,2,2);
plot(n1,x1);
xlabel("time");
ylabel("Amplitude");
grid;
title("Under sampling");
n2=-5:1:5;
x2=cos(2*pi*fm/fs2*n2);
subplot(2,2,3);
stem(n2,x2);
xlabel("Number of samples");
ylabel("Amplitude");
hold on;
subplot(2,2,3);
plot(n2,x2);
xlabel("time");
ylabel("Amplitude");
grid;
title("Uniform sampling");
n3=-20:1:20;
x3=cos(2*pi*fm/fs3*n3);
subplot(2,2,4);
stem(n3,x3);
xlabel("Number of samples");
ylabel("Amplitude");
hold on;
subplot(2,2,4);
plot(n3,x3);
xlabel("time");
ylabel("Amplitude");
grid;
title("Over sampling");



exp---999
butter BUTTER
clc
clear all
fs=100                  
f=5                     
t=5                     
n=[0:1/fs:t]            
x=2*sin(2*pi*f*n)       
subplot(3,1,1)
plot(n,x)
grid on
title('Sinusoidal signal');
z=awgn(x,1)          
subplot(3,1,2)
plot(n,z)
title('Sinusoidal signal with noise added');
o=1                
wc=2*pi*3.5/fs          
[b,a]=butter(1,wc,'low') 
iir=filter(b,a,z)
subplot(3,1,3)
plot(n,iir)
fvtool(b,a);   
%High - pass Filter
clc
clear all
fs=100                  
f=5                    
t=5                     
n=[0:1/fs:t]            
x=2*sin(2*pi*f*n)      
subplot(3,1,1)
plot(n,x)
grid on
title('Sinusoidal signal');
z=awgn(x,1)          
subplot(3,1,2)
plot(n,z)
title('Sinusoidal signal with noise added');
o=1              
wc=2*pi*3.5/fs           
[b,a]=butter(1,wc,'high') 
iir=filter(b,a,z)
subplot(3,1,3)
plot(n,iir)
fvtool(b,a);   



1010101010101010
PROCEDURE:
i) Analog filter: Band pass filter
clear all;

rp = input(‘pass ripple freq’);
rs = input(‘stop ripple freq’);
fp = input(‘pass band freq’);
fs = input(‘stop band freq’);
f = input(‘sample freq’);
w1 = 2*fp/f;
w2= 2*fs/f;
[nJ= buttord(w1 ,w2,rp,rs); wn= [w1,w2];
[b,a]= butter(n,wn,’bandpass’); w=0:0.1:pi;
[h,p]= freqz(b,a,w);
g= 20*log10(abs(h)); A=angle(h);
subplot (2,2,1); plot(p/pi,g); ylabel(‘amp’);
xlabel(‘ferq’);
title(‘amp,freq’);
subplot (2,2,2); plot(p/pi,A); xlabel(‘normal. freq’); ylabel(‘phase’);
title(‘normal. freq,phas&);



ii) IIR filter using impulse invariance: a sixth-order analog Butterworth lowpass filter to a
digital filter using impulse invariance
f = 2;
fs = 10;
[b,a] = butter(6,2*pi*f,'s');
[bz,az] = impinvar(b,a,fs);
freqz(bz,az,1024,fs)



iii) FIR filter: low pass filter
clc;
clear all;
wc=input('enter the value of cut off frequency');
N=input('enter the value of filter');
alpha=(N-1)/2;
eps=0.001;
%Rectangular Window
n=0:1:N-1;
hd=sin(wc*(n-alpha+eps))./(pi*(n-alpha+eps));
hn=hd
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h));
hold on
%Hamming Window
n=0:1:N-1;
wh=0.54-0.46*cos((2*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h),'ms');
hold off;
hold on
%Hanning Window
n=0:1:N-1;
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h),'blue');
hold off;
hold on
%Blackman Window
n=0:1:N-1; wh=0.42-0.5*cos((2*pi*n)/(N-1))+0.08*cos((4*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h),'green');
hold off;


iv) FIR filter: highpass filter
clc;
clear all;
wc=input('enter the value of cut off frequency');
N=input('enter the value of filter');
alpha=(N-1)/2;
eps=0.001;
%Rectangular Window
n=0:1:N-1;
hd=(sin(pi*(n-alpha+eps))-sin((n-alpha+eps)*wc))./(pi*(n-alpha+eps)); hn=hd
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h));
hold on
%Hamming Window
n=0:1:N-1;
wh=0.54-0.46*cos((2*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h),'ms');
hold off;
hold on
%Hanning Window
n=0:1:N-1;
wh=0.5-0.5*cos((2*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h),'blue');
hold off;
hold on
%Blackman Window
n=0:1:N-1; wh=0.42-0.5*cos((2*pi*n)/(N-1))-0.08*cos((4*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(w/pi,abs(h),'green');
hold off;



v)FIRfilter:bandstopfilter
clc;
Wc1=input('enter the value of Wc1=');
Wc2=input('enter the value of Wc2=');
N=input('enter the value of N=');
alpha=(N-1)/2;
eps=0.001;
%Rectangular Window
n=0:1:N-1; hd=(sin(Wc1*(n-alpha+eps))-sin(Wc2*(n-alpha+eps)*pi))./((nalpha+eps)*pi); hn=hd
W=0:0.01:pi;
h=freqz(hn,1,W);
plot(W/pi,abs(h));
hold on;
%Hamming Window
n=0:1:N-1;
Wn=0.54-0.46*cos((2*pi*n)/(N-1));
hn=hd.*Wn
W=0:0.01:pi;
h=freqz(hn,1,W);
plot(W/pi,abs(h),'green');
hold on;
%Hanning Window
n=0:1:N-1;
Wn=0.5-0.5*cos((2*pi*n)/(N-1));
hn=hd.*Wn
W=0:0.01:pi;
h=freqz(hn,1,W);
plot(W/pi,abs(h),'red');
hold off;
%Blackman Window
n=0:1:N-1;
wh=0.42-0.5*cos((2*pi*n)/(N-1))-0.08*cos((4*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(W/pi,abs(h),'green');
hold off;




vi) FIR filter: bandpass filter
clc;
Wc1=input('enter the value of Wc1=');
Wc2=input('enter the value of Wc2=');
N=input('enter the value of N=');
alpha=(N-1)/2;
eps=0.001;
%Rectangular Window
n=0:1:N-1;
 hd=(sin(Wc1*(n-alpha+eps))-sin(Wc2*(n-alpha+eps)*pi))./((n-alpha+eps)*pi);
hn=hd
W=0:0.01:pi;
h=freqz(hn,1,W);
plot(W/pi,abs(h));
hold on;
%HammingWindow
Wn=0.54-0.46*cos((2*pi*n)/(N-1));
hn=hd.*Wn
W=0:0.01:pi;
h=freqz(hn,1,W);
plot(W/pi,abs(h),'green');
hold on;
%Hanning Window
n=0:1:N-1;
Wn=0.5-0.5*cos((2*pi*n)/(N-1));
hn=hd.*Wn
W=0:0.01:pi;
h=freqz(hn,1,W);
plot(W/pi,abs(h),'red');
hold off;
%Blackman Window
n=0:1:N-1;
wh=042-0.5*cos((2*pi*n)/(N-1))-0.08*cos((4*pi*n)/(N-1));
hn=hd.*wh
w=0:0.01:pi;
h=freqz(hn,1,w);
plot(W/pi,abs(h),'green');
hold off;

