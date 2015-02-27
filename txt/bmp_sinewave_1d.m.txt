function bmp_sinewave_1d
%creates a graph of a 1D sinewave

Fs = 8000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.25;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
% Sine wave:
Fc1 = 15;                     % hertz
A1 = 1; %amplitude
Fc2 = 41;
A2 = 0;
Fc3 = 358;
A3 = 0;
x =  (A1* cos(2*pi*Fc1*t)) +  (A2* cos(2*pi*Fc2*t))+ (A3*   cos(2*pi*Fc3*t));
mx = max(x(:));
mn = min(x(:));
x = x  / (mx-mn);
% Plot the signal versus time:
figure;
set(gcf, 'color', [1 1 1])
plothandle = plot(t,x );
set( plothandle, 'Color', [ 0, 0, 0 ] );  % plot color to black
xlabel('time (in seconds)');
title('Signal versus Time');
axis([0,StopTime,-1.2,1.2])
axis off
set(gcf, 'Position', [30 30 800 300])
set(gca,'pos', [0 0 1 1]); 

 
