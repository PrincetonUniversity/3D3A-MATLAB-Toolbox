start3D3A

%% Create and Initialize
portno = 9000;
hudpr = dsp.UDPReceiver('LocalIPPort',portno);

%% Stream processing loop
fprintf('Listening on port %d...\n',portno);
tic;
while toc < 10
    % Read buffer from UDP
    udpBuff = step(hudpr);
    if(~isempty(udpBuff))
        lastBuff = udpBuff;
        disp('Message received.')
        [msg, dat] = readOSCmessage(udpBuff);
    end
end
disp('Done listening.')

%% Terminate
release(hudpr);