function [sig] = pcmread(url)
fileId = fopen(url,'r');
sig = fread(fileId,inf,'int16');
fclose(fileId);