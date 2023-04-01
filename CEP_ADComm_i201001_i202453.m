clc;close;clear;                                                      
%-------Creating the Audio-----------
%using sameple code given by instructor
Fs = 48000; %Sampling Frequency
recObj = audiorecorder; 
nBits = 16;
nChannels = 1;
ID = -1;
recObj = audiorecorder(Fs,nBits,nChannels,ID);

%message1
continue_enter = input('[Enter] Record Message 1'); %these inputs to break code in between
disp('Start speaking.');
recordblocking(recObj,5);
disp('End of Recording.');
play(recObj);
mySpeech = getaudiodata(recObj);
audiowrite('test1.wav',mySpeech,Fs);

%message2
continue_enter = input('[Enter] Record Message 2');           
disp('Start speaking.');                                    
recordblocking(recObj,5);                                  
disp('End of Recording.');
play(recObj);
mySpeech = getaudiodata(recObj);
audiowrite('test2.wav',mySpeech,Fs);

%%                                                          
%----------Read Files-------------

%read messages
figure (1)  
s1 = audioread('test1.wav');                                                                        
s2 = audioread('test2.wav');

%Play Message1
continue_enter = input('[Enter] Play Audio 1 ');
player6 = audioplayer(s1,44100);
playblocking(player6);
Am_1 = max(s1);
      
%Play Mesasge2
player7 = audioplayer(s2,44100);
continue_enter = input('[Enter] Play Audio 2 ');
playblocking(player7);
Am_2 = max(s2);


subplot(2,1,1)  %plot time domain of message signals                                            
plot(s1);                                                   
title('Message Signal 1');                                  
xlabel('t');                                                
ylabel('Amplitude');
subplot(2,1,2)                                                 
plot(s2);
title('Message Signal 2');
xlabel('t');
ylabel('Amplitude');
%%

%--------Parameters------------
guard = 1000;       %guard band
SNR = 2000;       %signal to noise ratio for awgn noise
%Carrier Frequencies for SSB
carrier1 = 6000;
carrier2 = 10000;
bandlimit_freq = 3000;
%Sampling Variables 
N = 5*Fs;   %number of samples    
n  = 0:N-1;	%number of samples vector    
Ts = 1/Fs;	%time period of signal
t  = (n*Ts)';   
%%    
%----------Transmitter---------
figure (2)                    
subplot(211);   %plot time domain and freq domain of message 1
plot(t,s1,'linewidth',2); 
title('Message Signal 1');
grid on;     
xlabel('t','fontweight','bold');
ylabel('message 1','fontweight','bold');
Xk2 = fft(s1)/N;                         
Y2 = abs(Xk2);                           
X2 = fftshift(abs(Xk2));                 
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));
subplot(212)                                   
plot(fx2,X2);                                  
grid on;                                                       
xlabel('-Fs/2 to Fs/2','fontweight','bold');                
ylabel('M1(f)','fontweight','bold');                        
                                                            
figure (3)                                                  
subplot(211);   %plot time domain and freq domain of message 2                                               
plot(t,s2,'linewidth',2);                                   
title('Message Signal 2');                                  
grid on;                                                        
xlabel('t','fontweight','bold');                               
ylabel('message 2','fontweight','bold');                       
Xk2 = fft(s2)/N;                                                  
Y2 = abs(Xk2);                                              
X2 = fftshift(abs(Xk2));                                    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));                 
subplot(212)                                                
plot(fx2,X2);                                               
grid on;                                                        
xlabel('-Fs/2 to Fs/2','fontweight','bold');                
ylabel('M2(f)','fontweight','bold');                        

%band limit signals
continue_enter = input('[Enter] Lowpass to band limit');                      
s1 = lowpass(s1, bandlimit_freq, Fs);         
s2 = lowpass(s2, bandlimit_freq, Fs);                                      

%-----------USSB Modulation----------

continue_enter = input('[Enter] SSB Modulation');
%message1        
s1_h = imag(hilbert(s1));
s1mod  = (1/2).*(s1.*cos(2*pi*carrier1*t) - s1_h.*sin(2*pi*carrier1*t));
s1mod = bandpass(s1mod, [6500 8500], Fs);

%mesasge2
s2_h = imag(hilbert(s2));
s2mod  = (1/2).*(s2.*cos(2*pi*carrier2*t) - s2_h.*sin(2*pi*carrier2*t));
s2mod = bandpass(s2mod, [10500 13000], Fs);

%Plot Modulated M1
figure (4)                                                  
subplot(211);                                               
Xk2 = fft(s1mod)/N;                                             
Y2 = abs(Xk2);                                              
X2 = fftshift(abs(Xk2));                                    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));             
plot(fx2,X2);                                               
title('USSB Modulated Signals');                            
xlabel('-Fs/2 to Fs/2','fontweight','bold');              

%Plot Modulated M2
subplot(212)                                                
Xk2 = fft(s2mod)/N;                                             
Y2 = abs(Xk2);                                              
X2 = fftshift(abs(Xk2));                                    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));             
plot(fx2,X2);                                               
grid on;                                                    
xlabel('-Fs/2 to Fs/2','fontweight','bold');                                                                          

%------------Multiplexing------------

continue_enter = input('[Enter] Multiplexing The Signals'); 
FDM = s1mod+s2mod;

figure (5)                                                   
subplot(211)    %plotting multiplexed signal                                              
plot(t,FDM,'linewidth',2);                                  
title('Multiplexed Signal');                                
grid on;                                                    
xlabel('t','fontweight','bold');                                
ylabel('FDM','fontweight','bold');                          
                                                            
subplot(212)                                                       
Xk2 = fft(FDM)/N;                                               
Y2 = abs(Xk2);                                              
X2 = fftshift(abs(Xk2));                                    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));             
plot(fx2,X2);                                               
title('Frequency Domain');                                  
grid on;                                                       
xlabel('-Fs/2 to Fs/2','fontweight','bold');                                                                   
continue_enter = input('[Enter] Frequency Modulation & Deomodulation');       
                                                            
%--------Freq Modulation---------

%FM Carrier Variables
FM_Fc = 1000000;    %carrier wave frequency
FM_Fs = 3*FM_Fc;    %sampling rate 
FM_BW = 200e3;  %bandwidth
freqdev = FM_BW/2;  %as it is wideband fm -> BW = 2freqdev                       

%FM Modulation
y = fmmod(FDM,FM_Fc,FM_Fs,freqdev); %using freq modulation function
  
figure (6)             
subplot(211)
plot(t,y,'linewidth',2);
title('Frequency Modulated Signal');
grid on;                            
xlabel('t','fontweight','bold');    
ylabel('FM_Modulation','fontweight','bold');
                                            
subplot(212)                                
Xk2 = fft(FDM)/N;
Y2 = abs(y);     
X2 = fftshift(abs(Xk2));    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2)); 
plot(fx2,X2);                
title('Frequency Domain');                      
grid on;                                        
xlabel('-Fs/2 to Fs/2','fontweight','bold');    

%noise
y = awgn(y, SNR); 

figure (7)                                                  
plot(t,y,'linewidth',2);                                    
title('Adding Noise in Wave');                              
grid on;                                                    
xlabel('t','fontweight','bold');                            
ylabel('FM_Modulation + Noise','fontweight','bold');        

%FM Demodulation
FDM = fmdemod(y,FM_Fc,FM_Fs,freqdev);   %demodulating the FM signal
                                    
figure (8)                          
subplot(211)                        
plot(t,FDM,'linewidth',2);          
title('Frequency Demodulated Signal');
grid on;                              
xlabel('t','fontweight','bold');      
ylabel('FM_Demodulation','fontweight','bold');              
                                                                
subplot(212)                                                
Xk2 = fft(FDM)/N;                                           
Y2 = abs(FDM);                                             
X2 = fftshift(abs(Xk2));                                                                            
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));             
plot(fx2,X2);                                                   
title('Frequency Domain');                                 
grid on;                                                   
xlabel('-Fs/2 to Fs/2','fontweight','bold');                


%Demultiplexing Signals using bandpass filters
continue_enter = input('[Enter] Demultiplexing by filtering bands');                                                                                 
demuxs1 = bandpass(FDM, [6500 8500], Fs);                                
demuxs2 = bandpass(FDM, [10500 13000], Fs);

figure (9)                                                  
                                                            
subplot(221)                                                
plot(t,demuxs1,'linewidth',2);                              
title('Demultiplexed M1');                                  
grid on;                                                        
xlabel('t','fontweight','bold');                                
ylabel('FDM','fontweight','bold');                          
                                                            
subplot(222)                                                
Xk2 = fft(demuxs1)/N;                                       
Y2 = abs(Xk2);                                              
X2 = fftshift(abs(Xk2));                                    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));             
plot(fx2,X2);                                               
title('Frequency Domain');                                  
grid on;                                                    
xlabel('-Fs/2 to Fs/2','fontweight','bold');                
                                                            
subplot(223)                                                
plot(t,demuxs2,'linewidth',2);                              
title('Demultiplexed M2');                                  
grid on;                                                    
xlabel('t','fontweight','bold');                                
ylabel('FDM','fontweight','bold');                          
                                                            
subplot(224)                                                
Xk2 = fft(demuxs2)/N;                                       
Y2 = abs(Xk2);                                              
X2 = fftshift(abs(Xk2));                                    
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));                 
plot(fx2,X2);                                               
title('Frequency Domain');                                  
grid on;                                                       
xlabel('-Fs/2 to Fs/2','fontweight','bold');                    
                                              
%USSB Demodulation

%message1
carrier1_demod = cos(2*pi*carrier1*t);
demods1 = demuxs1 .* carrier1_demod;
demods1 = lowpass(demods1, 3000, Fs);                       

carrier2_demod = cos(2*pi*carrier2*t);
demods2 = demuxs2 .* carrier2_demod;                                                                
demods2 = lowpass(demods2, 3000, Fs);                       
                                                            
continue_enter = input('[Enter] Recovered Message 1');                                                                            
player4 = audioplayer(demods1,44100);                       
playblocking(player4);                                          
                                                            
continue_enter = input('[Enter] Recovered Message 2');                                                                 
player5 = audioplayer(demods2,44100);                       
playblocking(player5);                                      
                                                            
figure (10)                                                 
subplot(2,1,1)                                              
plot(demods1);                                              
title('Recovered Message Signal 1');                        
xlabel('t');                                                
ylabel('Amplitude');                                        
subplot(2,1,2)                                              
plot(demods2);                                              
title('Recovered Message Signal 2');                        
xlabel('t');                                                
ylabel('Amplitude');                                        
%%----end of code------


