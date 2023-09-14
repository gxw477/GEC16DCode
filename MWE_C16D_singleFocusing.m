%% MWE 

clear all
% Specify system parameters

P.numRays = 128;  % no. of raylines to program
P.numTx = 128;
P.startDepth = 5; % startDepth and endDepth in wavelength
P.endDepth = 200;

P.TXFocus = [200 ];%:-50:100;
P.numFz = length(P.TXFocus);
P.wait = 200; % 50 msec pause

Resource.Parameters.numTransmit = 128; % no. of transmit channels
Resource.Parameters.numRcvChannels = 128; % change to 64 for Vantage 64 system
Resource.Parameters.connector = 1; % trans. connector to use (V 256).
Resource.Parameters.speedOfSound = 1540; % speed of sound in m/sec
Resource.Parameters.simulateMode = 1; % runs script in simulate mode
% Specify media points
%Media.MP(1,:) = [0,0,100,1.0]; % [x, y, z, reflectivity]
% Specify Trans structure array.
Trans.name = 'GEC1-6D';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C1-6D transducer is 'known' transducer so we can use computeTrans.


% Specify Resource buffers.
samplesPerWave = 2; 
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2^18;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 2; % minimum size is 1 frame.
% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1]; % A, B, C, D

%% Pdata  Tings

radius = Trans.radius;
dr = 1;%Radial increment

scanangle = Trans.numelements*Trans.spacing/radius;
dtheta = scanangle/P.numRays;
theta = -(scanangle/2) + 0.5*dtheta; % angle to left edge from centerline
Angle = theta:dtheta:(-theta);

P.dtheta = dtheta;

m=1;

for n = 1:P.numRays
    PData.Region(n+m) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta,...
                       'steer',Angle(n),...
                       'andWithPrev',1));
end

%% Transmit structure

for iTransmit = 1:P.numRays
    PData(1).Region(iTransmit) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta,...
                       'steer',Angle(iTransmit)));
end


TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.numRays*P.numFz);

% Determine TX aperture based on focal point and desired f number.
txFNum = 3.5;  % set to desired f-number value for transmit (range: 1.0 - 20)
txNumEl = zeros(1,P.numFz);

for j = 1:P.numFz
    txNumEl(j)=round((P.TXFocus(j)/txFNum)/Trans.spacing/2); % no. of elements in 1/2 aperture.
    if txNumEl(j) > (Trans.numelements/2 - 1), txNumEl(j) = floor(Trans.numelements/2 - 1); end
end

for iTransmit = 1:P.numRays   % 128 transmit events
    
    %pick the starting point in the TX struct array. There are P.numFz per
    %lateral position
    k = P.numFz*(iTransmit-1);

    %cycle through focal depths (from deepest to shallowest)
    for n = 1:P.numFz
        % Set transmit Origins to positions of elements.
        TX(n+k).Origin = [radius*sin(Angle(iTransmit)), 0.0, radius*cos(Angle(iTransmit))-radius];
        % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.

        ce = round(1+Trans.numelements*(Angle(iTransmit) - theta)/(-2*theta));
        % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
        lft = round(ce - txNumEl(n)/2);
        if lft < 1, lft = 1; end;
        rt = round(ce + txNumEl(n)/2);
        if rt > Trans.numelements, rt = Trans.numelements; end;
        TX(n+k).Apod(lft:rt) = 1.0;
        TX(n+k).focus = P.TXFocus(n);
        TX(n+k).Delay = computeTXDelays(TX(n+k));
    
    end
end




%% Receive structure 
% Specify TGC Waveform structure.
TGC(1).CntrlPts = [500,590,650,710,770,830,890,950];
TGC(1).rangeMax = 200;
TGC(1).Waveform = computeTGCWaveform(TGC);


maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
RcvZone = maxAcqLength;

%theta = -63.5*Trans.spacing/radius; % angle to element 1 of transducer
thetaLeft = -(Trans.numelements/2 - round(P.numTx/8))*Trans.spacing/radius; % angle to lft edge of scan

PData.PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData.PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData.PDelta(1));
PData.Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData.Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData.Origin(1,2) = 0;
PData.Origin(1,3) = ceil(radius * cos(thetaLeft)) - radius - 5;


m = m + P.numRays;
for n = 1:P.numRays
    if n<10
        steer = -((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        steer = -((P.numRays-n)/10)*P.dtheta;
    else
        steer = -P.dtheta;
    end
    d = radius*tan(-steer);
    b = sqrt(radius*radius + d*d);
    c = sqrt(2*radius*radius + d*d - 2*radius*b*cos(-steer));
    PData.Region(n+m) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[c,0,-radius],...
                       'r1',b+P.startDepth,...
                       'r2',b+P.endDepth,...
                       'angle',P.dtheta*6,...
                       'steer',(Angle(n)+steer),...
                       'andWithPrev',1));
end

m = m + P.numRays;
for n = 1:P.numRays
    if n<=10
        steer = ((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        steer = ((P.numRays-n)/10)*P.dtheta;
    else
        steer = P.dtheta;
    end
    d = radius*tan(steer);
    b = sqrt(radius*radius + d*d);
    c = sqrt(2*radius*radius + d*d - 2*radius*b*cos(steer));
    PData.Region(n+m) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[-c,0,-radius],...
                       'r1',b+P.startDepth,...
                       'r2',b+P.endDepth,...
                       'angle',P.dtheta*6,...
                       'steer',(Angle(n)+steer),...
                       'andWithPrev',1));
end
   
PData.Region = computeRegions(PData);

%% Recon struct

Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, P.numRays);

% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = 1+P.numFz*(j-1); % only need the 1st tx out of each 3 tx (for different focal zones) as they have the same origin and characteristics (focus)
    ReconInfo(j).rcvnum = 1+P.numFz*(j-1);
    ReconInfo(j).regionnum = j;
end



%% Receive struct 

Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'bufnum', 1, ...
                        'TGC',1, ....
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'BS100BW',...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1,P.numFz*P.numRays*Resource.RcvBuffer(1).numFrames);

m = P.numRays;
for i = 1:Resource.RcvBuffer(1).numFrames
    
    k = P.numRays*P.numFz*(i-1);
    Receive(k+1).callMediaFunc = 1;
    
    for j = 1:P.numRays
        
        w =  P.numFz*(j-1);

        for z = 1:P.numFz
            %Set max acq length for first TX in each set of focal zones
            if z == 1
                Receive(k+w+z).endDepth = maxAcqLength;
            else
                Receive(k+w+z).endDepth = RcvZone(z);
            end
            Receive(k+w+z).Apod(:) = 1.0;
            Receive(k+w+z).framenum = i;
            Receive(k+w+z).acqNum = j;      % two acquisitions per frame
         
        end
    end

end


Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 10;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 10;
Resource.DisplayWindow(1).Title = 'GE4CDWideBeam';
Resource.DisplayWindow(1).pdelta = 0.5;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 10;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

%% External processing 

% Specify an external processing event.

pers = 0;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',0,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',10,...            % pgain is image processing gain
    'persistMethod','simple',...
    'persistLevel',pers,...
    'reject',2,...      % reject level
    'interpMethod','4pt',...  %method of interp. (1=4pt)
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','full',... 
    'display',1,...      % display image after processing
    'displayWindow',1};


points_noSpeckle
Media.function = 'movePoints';

%% Events
% Specify SeqControl structure arrays.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;  % 200 usec between ray lines; should be > 2*Receive.enDepth*T(period)
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(80000 - P.numRays*P.numFz*200); % 12.5 frames per second, fz=3
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'jump'; % Jump back to start.
SeqControl(4).argument = 1;

iJumps = 1;
iRecons = 1;
nEventsTotal =  Resource.RcvBuffer(1).numFrames*P.numFz*P.numRays + ...
    Resource.RcvBuffer(1).numFrames*iRecons + iJumps;
            

Event = repmat(struct(...
                'info',[],...
                'tx',[],...
                'rcv', [],...
                'recon',[],...
                'process',[],...
                'seqControl',[]), ...
                1,nEventsTotal);

nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = P.numRays*P.numFz*(i-1);
    for j = 1:P.numRays                 % Acquire all ray lines for frame
        w = P.numFz*(j-1);
        for z = 1:P.numFz
            Event(n).info = 'Acquire ray line';
            Event(n).tx = z+w;
            Event(n).rcv = k+z+w;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 1;
            n = n+1;
        end
    end
    % Replace last event's SeqControl with inter-frame timeToNextAcq and transfer to host.
    Event(n-1).seqControl = [2,nsc];
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    if floor(i/4) == i/4     % Exit to Matlab every 4th frame reconstructed
        Event(n).seqControl = 3;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;


% Save all the structures to a .mat file.



UI(1).Control = {'UserB1','Style','VsPushButton','Label','saveRF'};
UI(1).Callback = text2cell('%-UI#1Callback');


save('MWE/SetupData/sFocus');

return

%-UI#1Callback - save RF
if evalin('base','freeze')==0   % no action if not in freeze
    msgbox('Please freeze VSX');
    return
end

Control.Command = 'copyBuffers';
runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.

RFfilename = ['RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];

RcvLastFrame = size(RcvData,3);
if (~evalin('base','simButton'))
    RcvLastFrame = Resource.RcvBuffer(1).lastFrame;
end

[fn,pn,~] = uiputfile('*.mat','Save RF data as',RFfilename);
if ~isequal(fn,0) % fn will be zero if user hits cancel
    fn = strrep(fullfile(pn,fn), '''', '''''');
    save(fn,'RcvData','RcvLastFrame','-v7.3');
    fprintf('The RF data has been saved at %s \n',fn);
else
    disp('The RF data is not saved.');
end

return
%-UI#1Callback




