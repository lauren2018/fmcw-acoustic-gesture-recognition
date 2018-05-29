clear
Fs = 48000;
socket = tcpip('0.0.0.0',12345,'NetworkRole', 'server');
set(socket,'InputBufferSize', Fs*10)
fopen(socket);
'connected'
dataAll = [];
[b, a] = butter(10, 15500/(Fs/2), 'high');

Fs = 48000;
initFreq = 16e3;
finlFreq = 22e3;	
timeInterv = 0.02;
t = 0:1/Fs:timeInterv;
chirpSlice = chirp(t,initFreq,timeInterv,finlFreq); 
N = length(chirpSlice) ;
H = hanning(N);
hanChirp = chirpSlice.*H';
finChirp=[];
for i=1:248
    finChirp=[finChirp,hanChirp];
end
chirp10x = finChirp;
syn = 0;
offset = 0;

while(1)
    while(socket.BytesAvailable==0)
    end
    socket.BytesAvailable
    data = fread(socket,socket.BytesAvailable);
    data = data';
    dataint16 = typecast(uint8(data),'int16');
    dataAll = [dataAll dataint16];
    if length(dataAll) < 5*Fs
        continue
    end
    data = double(dataAll(end-5*Fs+1 : end));
%     data = double(dataAll(end-5*Fs : end));
    data = filtfilt(b, a, data);
    if (~syn)
        [f, l] = xcorr(data, chirp10x);
        [~, I] = max(f);
        offset = l(I);
        if(offset <= 0 || offset+length(chirp10x)+1 > length(data))
    %     if(true)
            plot(data)
            drawnow
%             pause(0.5)
            continue
        else
            offset = l(I); 
            syn = 1;
        end
    end
    rcvChirp = double(data(offset+1 : offset+length(chirp10x)));
    deChirp = rcvChirp .* chirp10x;
    sig = {};
    for i = 1 : length(hanChirp) : length(chirp10x)
       f = abs(fft(deChirp(i:i+length(hanChirp)-1), 8192*2)); 
%        plot(f(1:1000))
%        drawnow
       sig = [sig, f(1:1000)];
    end
    
    idx = [];
    for i = 1 : length(sig)-1
        t = sig{i+1}-sig{i};
        t = t - mean(t);
        t = abs(t);
        [pks, locs] = findpeaks(t);
        [M, l] = max(pks);
%         [M, j] = max(t);
        idx = [idx, locs(l)];
    end
%     idx = idx-mean(idx);
%     plot(idx)
    plot(t)
    drawnow

    win = 50;
    dist = [];
    for i = 1:length(idx)-win
        dist = [dist, sum(idx(i:i+win-1))];
    end
%     plot(dist)
%     drawnow
%     ylim([-1000, 1000])
%     pause(0.5)
end
