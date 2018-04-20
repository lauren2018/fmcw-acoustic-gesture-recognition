% Author: Wentao Xie, Meng Zhou and Xiaotong Zhang
clear;
close all;

%% --------------- Constants & configuration ------------------

Fs = 48000;
initFreq = 17e3;
finlFreq = 19.5e3;	
timeInterv = 0.04;

% -------------------------------------------------------------
%% --------------- A frame of chirp signal --------------------

t = 0:1/Fs:timeInterv;
chirpSlice = chirp(t,initFreq,timeInterv,finlFreq); 

% -------------------------------------------------------------
%% ----------------- Series of chirp sigs ---------------------

N = length(chirpSlice) ;
H = hanning(N);
hanChirp = chirpSlice.*H';
finChirp=[];
for i=1:250
    finChirp=[finChirp,hanChirp];
end

% -------------------------------------------------------------
%% ---------------------- Write / read ------------------------

douTrackSig = [finChirp', zeros(1, length(finChirp))'];
% audiowrite('test.wav', douTrackSig, Fs);
rcvChirp = pcmread('move2.pcm');

% -------------------------------------------------------------
%% ----------------------- Dechirp ----------------------------

[f,idx] = xcorr(rcvChirp, finChirp);
[~,I] = max(f);
offset = idx(I);
rcvChirp = rcvChirp(offset+1 : offset+length(finChirp));
deChirp = finChirp .* rcvChirp';
for i = 1:length(hanChirp):length(finChirp)
      f = abs(fft(deChirp(i:i+length(hanChirp)-1)));
      plot(f(1:100))
      ylim([0, 290000])
      xlim([0, 100])
      drawnow
end

% -------------------------------------------------------------