clc
clear
t=0:0.01:1;
a=2;
b=a*sin(2*pi*2*t);
subplot(3,3,1);
stem(t,b);
xlabel('Time');
ylabel('Amplitude');
title('Sine wave');
t=0:0.01:1;
s=2;
b=s*cos(2*pi*2*t);
subplot(3,3,2);
stem(t,b);
xlabel('Time');
ylabel('Amplitude');
title('Cosine wave');
n=
-5:5;
a=[zeros(1,5),ones(1,1),zeros(1,5)];
subplot(3,3,3);
stem(n,a);
xlabel('Time');
ylabel('Amplitude');
title('Unit impulse');
3
n=
-5:5;
a=[zeros(1,5),ones(1,6)];
subplot(3,3,7);
stem(n,a);
xlabel('Time');
ylabel('Amplitude');
title('Unit step');
t=0:0.01:1
;
a=2;
b=a*square(2*pi*2*t);
subplot(3,3,8);
stem(t,b);
xlabel('Time');
ylabel('Amplitude');
title('Square wave');
t=0:0.01:1;
r=2;
b=r*exp(2*pi*2*t);
subplot(3,3,9);
stem(t,b);
xlabel('Time');
ylabel('Amplitude');
title('Exponential Wave');






DFT DFT
clc;
clear all;
close all;
x=input('Enter input sequence:');
N=length(x);
disp(‘Length:);
disp(N)
subplot(3,2,1);
stem(x);
xlabel('n value');
ylabel('Amp');
title('Input value');
y=zeros(1,N);
for k=0:N-1
 for n=0:N-1
 y(k+1)=y(k+1)+x(n+1)*exp((-2*i*pi*k*n)/N);
 end
end
disp(‘DFT:’);
disp('y=');
disp(y)
subplot(3,2,2);
stem(y);
xlabel('n value');
ylabel('Amp');
title('DFT value');
magnitude=abs(y);
7
subplot(3,2,3);
stem(magnitude);
xlabel('n value');
ylabel('Amp');
title('Magnitude');
z=phase(y);
subplot(3,2,4);
stem(z);
xlabel('n value');
ylabel('Amp');
title('phase');
%IDFT
N=length(y);
m=zeros(1,N);
for n=0:N-1
 for k=0:N-1
 m(n+1)=m(n+1)+((1/N)*(y(k+1)*exp((i*2*pi*k*n)/N)));
 end
end
disp(IDFT:);
disp('m=')
disp(m);
subplot(3,2,5);
stem(m);
xlabel('n value');
ylabel('Amp');
title('IDFT');


exp 3 linear convulation:
Clc;
clear all;
close all;
x1=input('Enter the first sequence:');
10
subplot(3,1,1);
stem(x1);
xlabel('Time');
ylabel('Amplitude');
title('Plot the first sequence');
x2=input('Enter second sequence:');
subplot(3,1,2);
stem(x2);
xlabel('Time');
ylabel('Amplitude');
title('Plot of second sequence');
f=conv(x1,x2);
disp('Output of linear convolution is');
disp(f); %disp is used for display the message in command window
subplot(3,1,3);
stem(f);
xlabel('Time index n');
ylabel('Amplitude f');
title('Linear convolution of sequence');




Auto correlation

clc;
close all;
clear all;
%input sequences
x=input('Enter input sequence:')
subplot(1,2,1);
stem(x);
13
xlabel('n');
ylabel('x(n)');
title('Input Sequence');
% auto correlation of input sequence
z=xcorr(x,x);
disp('The values of z are = ');
disp(z);
subplot(1,2,2);
stem(z);
xlabel('n');
ylabel('z(n)');
title(' Auto correlation of input sequence'); 






FFFFFFFTTTTTTTTTTT
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
title(" Phase figure")






Interpolationnnnnnnnnnnn
clc;
clear all;
close all;
F=input('Enter the frequency of the signal: ');
P=input('Enter the interpulator factor: ');
18
N=input('Enter the length of the signal: ');
t=0:1:N-1;
X = sin(2*3.14*F*t);
subplot(2,1,1);
stem(X);
title("Original signal");
i=interp(X,P);
subplot(2,1,2);
stem(i);
title("Interpolated signal");












NNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDFFFFFFFFFFFtttttttttttttttt
close all;
clear all;
clc;
x1 = input('Enter the sequence: ');
n = input('Enter the length ');
20
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
%phase plot
an = angle(fft(x1,n));
disp(an);
subplot(2,2,2);
stem(N,an);
xlabel('Length');
ylabel('Phase of x(k)');
title('phase spectrun');



sampling
clc;
clear all;
close all;
t=-100:01:100;
fm=0.02;
x=cos(2*pi*t*fm);
subplot(2,2,1);
23
plot(t,x);
xlabel('time in sec');
ylabel('x(t)');
title('Continuous Time Signal');
fs1=0.02;
n=-2:2;
x1=cos(2*pi*fm*n/fs1);
subplot(2,2,2);
stem(n,x1);
hold on
subplot(2,2,2);
plot(n,x1,':');
title('Discrete Time Signal x(n) with fs<2fm');
xlabel('n');
ylabel('x(n)');
fs2=0.04;
n1=-4:4;
x2=cos(2*pi*fm*n1/fs2);
subplot(2,2,3);
stem(n1,x2);
hold on
subplot(2,2,3);
plot(n1,x2,':');
title('Discrete Time Signal x(n) with fs=2fm');
xlabel('n');
ylabel('x(n)');
24
n2=-50:50;
fs3=0.5;
x3=cos(2*pi*fm*n2/fs3);
subplot(2,2,4);
stem(n2,x3);
hold on
subplot(2,2,4);
plot(n2,x3,':');
xlabel('n');
ylabel('x(n)');
title('Discrete Time Signal x(n) with fs>2fm');











Exp-99999999
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
27
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
