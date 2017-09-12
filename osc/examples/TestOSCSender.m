start3D3A

%% Create and Initialize
portno = 9000;
hudps = dsp.UDPSender('RemoteIPPort',portno);
fprintf('Port set to %d.\n',portno);

%% Write message to UDP
msg = '/1/label1';
dat = {int32(4),4,'test'};
udpBuff = writeOSCmessage(msg, dat);
step(hudps,udpBuff);
disp('Message sent.')

%% Terminate
release(hudps);