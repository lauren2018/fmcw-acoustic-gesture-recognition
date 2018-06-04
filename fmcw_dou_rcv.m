% Author: Wentao Xie, Meng Zhou and Xiaotong Zhang
N = 8192*2;
res = (1*Fs) / N;

%% ------------------- Data pre-processing ----------------------

rcvSig = pcmread("e6.pcm")';
[b, a] = butter(10, 16000/(Fs/2), 'high');
rcvFiltered = filtfilt(b, a, rcvSig);

[f, l] = xcorr(rcvFiltered, model);
[~, I] = max(f);
offset = l(I);
rcvChirp = rcvFiltered(offset+1+1 : offset+length(model)+1);

rcvLeft = []; rcvRight = [];
for i = 1 : 1922 : length(rcvChirp)
    rcvLeft = [rcvLeft, rcvChirp(i : i+961-1)];
end
for i = 962 : 1922 : length(rcvChirp)
    rcvRight = [rcvRight, rcvChirp(i : i+961-1)];
end

%% ---------------------- Dechirp -------------------------

deChirpLeft = rcvLeft .* modelLeft;
deChirpRight = rcvRight .* modelRight;
Left = {}; Right = {};
for i = 1 : length(hanChirpLeft) : length(finChirpLeft)
    if i+length(hanChirpLeft)-1 > length(deChirpLeft)
        break
    end
    fLeft = abs(fft(deChirpLeft(i : i+length(hanChirpLeft)-1), N));
    fRight = abs(fft(deChirpRight(i : i+length(hanChirpRight)-1), N));
%     subplot(2, 1, 1)
%     plot(fLeft(1:1000))
%     subplot(2, 1, 2)
%     plot(fRight(1:1000))
%     drawnow
    Left = [Left, fLeft(1:1000)];
    Right = [Right, fRight(1:1000)];
end

%% ---------------------- Analysis -------------------------

idxLeft = []; idxRight = [];
for i = 1 : length(Left)-1
    dLeft = abs(Left{i+1} - Left{i});
    dRight = abs(Right{i+1} - Right{i});
    
%     subplot(2, 1, 1)
%     plot(dLeft)
%     subplot(2, 1, 2)
%     plot(dRight)
%     drawnow
%     pause

    if max(dLeft) < 8000 && i > 2
        idxLeft = [idxLeft, idxLeft(end)];
    else
        [pks, locs] = findpeaks(dLeft);
        [~, l] = max(pks);
        deltaF = locs(l)*res;
        deltaT = deltaF / k;
        idxLeft = [idxLeft, deltaT*vs];
        %     [~, jLeft] = max(dLeft);
        %     idxLeft = [idxLeft, jLeft/30000*340*100];
    end
    
    if max(dRight) < 8000 && i > 2
        idxRight = [idxRight, idxRight(end)];
    else
        [pks, locs] = findpeaks(dRight);
        [~, l] = max(pks);
        deltaF = locs(l)*res;
        deltaT = deltaF / k;
        idxRight = [idxRight, deltaT*vs];
        %     [~, jRight] = max(dRight); 
        %     idxRight = [idxRight, jRight/30000*340*100];
    end
end

% win = 50;
% distLeft = []; distRight = [];
% P = {};
% for i = 1 : length(idxLeft)-win
%     distLeft = [distLeft, mean(idxLeft(i : i+win-1))];
%     distRight = [distRight, mean(idxRight(i : i+win-1))];
%     P = [P, findPoint(distLeft(end), distRight(end))];
% end

distLeft = medfilt1(idxLeft, 10);
% disLeft = idxLeft;
distRight = medfilt1(idxRight, 10);
% disRight = idxRight;
P = {};
for i = 1 : length(idxLeft)
    [x, y] = findPoint4(distLeft(i), distRight(i));
    P = [P, [real(x), real(y)]];
end

figure
plot(idxLeft)
hold on
plot(idxRight)
figure
plot(distLeft)
hold on
plot(distRight)

figure
for i = 1 : length(P)
   plot(P{i}(1), P{i}(2), '.k')
   xlim([0, 0.009])
   ylim([-0.3, 0.3])
   hold on
%    drawnow
end
hold off

function [xp, yp] = findPoint(dl, dr)
opt = optimset('Display', 'off', 'Algorithm','Levenberg-Marquardt'); 
fun = @intEllipse;
x0 = [0, 0];
[xp, yp] = fsolve(fun, x0, opt);
    function F = intEllipse(X)
        x = X(1);
        y = X(2);
        F(1) = sqrt(x^2+y^2) + sqrt((x-3.5)^2+y^2) - dl;
        F(2) = sqrt((x-13.5)^2+y^2) + sqrt((x-3.5)^2+y^2) - dr;
        F(3) = (y+27/7*x-27/2)/abs(y+27/7*x-27/2) + 1;
    end
end

function [xp, yp] = findPoint2(dl, dr)
opt = optimset('Display', 'off', 'Algorithm','Levenberg-Marquardt'); 
fun = @intEllipse;
x0 = [1, 1];
[xp, yp] = fsolve(fun, x0, opt);
    function F = intEllipse(X)
        a = dl / 2;
        c = 13.5 / 2;
        b = sqrt(a^2 - c^2);
        r = dr / 2;
        x = X(1);
        y = X(2);
%         F(1) = y - sqrt(dr^2/4 - (x-6.75)^2);
%         F(2) = y - 1/2 * sqrt((dl^2-182.25)*(1-4*x^2/dl^2));
        F(1) = x^2/a^2 + y^2/b^2 - 1;
        F(2) = x^2 + y^2 - r^2;
        F(3) = y / abs(y) - 1;
    end
end

function [xp, yp] = findPoint3(dl, dr)
a = dl / 2;
c = 13.5 / 2;
b = sqrt(a^2 - c^2);
r = dr / 2;
syms x y
eqns = [x^2/a^2 + y^2/b^2 == 1, (x-c)^2+y^2 == r^2];
rs = solve(eqns, x, y);
xp = 1;
yp = 1;
end

function [xp, yp] = findPoint4(d1, d2)
L1 = 0.135;
L2 = 0;
xp = sqrt((d1^2-L1^2)*(d2^2-L2^2)*((L1+L2)^2-(d1-d2)^2)) / (2*d2*13.5);
yp = (d2*L1^2-d1*L2^2-d1^2*d2+d2^2*d1) / (2*(d1*L2+d2*L1));
end