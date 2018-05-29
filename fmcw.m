% Author: Wentao Xie, Meng Zhou and Xiaotong Zhang
clear;
% close all;

%% --------------- Constants & configuration ------------------

Fs = 48000;
initFreq = 17e3;
finlFreq = 23e3;	
timeInterv = 0.02;

% -------------------------------------------------------------
%% --------------- A frame of chirp signal --------------------

t = 0:1/Fs:timeInterv;
chirpSlice = chirp(t,initFreq,timeInterv,finlFreq); 

% -------------------------------------------------------------
%% ----------------- Series of chirp sigs ---------------------

N = length(chirpSlice) ;
H = hanning(N);
hanChirp = chirpSlice.*H';
% hanChirp = chirpSlice;
finChirp=[];
for i=1:1000
    finChirp=[finChirp,hanChirp];
end

% -------------------------------------------------------------
%% ---------------------- Write / read ------------------------

douTrackSig = [ zeros(1, length(finChirp))', finChirp' ];
% audiowrite('gesture.wav', douTrackSig, Fs);
rcvChirp = pcmread('t2.pcm');
[b, a] = butter(10, 16000/(Fs/2), 'high');
rcvChirp = filtfilt(b, a, rcvChirp);
% -------------------------------------------------------------
%% ----------------------- Dechirp ----------------------------
locs = [];
sig = {};
S = [];
[f,idx] = xcorr(rcvChirp, finChirp);
[~,I] = max(f);
offset = idx(I);
rcvChirp = rcvChirp(offset+1+1 : offset+length(finChirp)+1);
deChirp = finChirp .* rcvChirp';
% figure;
for i = 1:length(hanChirp):length(finChirp)
      f = abs(fft(deChirp(i:i+length(hanChirp)-1), 8192 *2));
%       plot(f(1:1000))
%       ylim([0, 250000])
%       xlim([0, 1000])
%       drawnow
%       pause(0.1)
      sig = [sig, f(1:1000)];
      S = [S, f(1:100)'];
      [pks, lcs] = findpeaks(f(1:1000));
%       locs = [locs, lcs(1)];
end

% -------------------------------------------------------------
idx = [];
D = [];
% figure;
for i = 1 : length(sig)-1
    t = sig{i+1}-sig{i};
    t = abs(t - mean(t));
    D = [D, t(1:20)'];
    win = 5;
%     plot(t);
%     ylim([0, 35000])
%     drawnow 
%     pause(0.2)
%     [~, j] = sort(t);
    [~, j] = max(t);
    [pks, locs] = findpeaks(t);
    [~, l] = max(pks);
    idx = [idx, locs(l)];
   
end
figure; plot(idx)
win = 50;
ii = [];
for i = 1:length(idx)-win
    ii = [ii, sum(idx(i:i+win-1))];
end
figure
plot(ii)
ii=table(ii);
