function [Stublinked, BlinkTimes] = Stublinks60(inputMat, SAMPLE_RATE)

%SAMPLE_RATE = 1000;
DURATION = size(inputMat,2)./SAMPLE_RATE; %12000/1000 = 12

%Stublinked = zeros(size(inputMat)); %size 2 x 12000
%BlinkTimes = zeros(size(inputMat)); %size 2 x 12000

%Stublinked = zeros(1,size(inputMat,2)+1); %size 2 x 12000
%BlinkTimes = zeros(1,size(inputMat,2)+1); %size 2 x 12000

% for each row (each eye) 
for i=1:size(inputMat,1)
    %resize
%     if mod(i,100) == 0
%         disp(i);
%     end
    
    % downsample to 720 (60Hz) using imresize 
    currentTrial = imresize(inputMat(i,:),[1 DURATION*60]); % currentTrial is size 1 x 720
    
    try
        output=stublinks(currentTrial);
        
        %resize output data
        outputData = imresize(output.NoBlinksUnsmoothed,[1 DURATION*SAMPLE_RATE]);
        %outputData = imresize(output.NoBlinks,[1 DURATION*SAMPLE_RATE]);
        outputTimes = imresize(output.BlinkTimes',[1 DURATION*SAMPLE_RATE]);
        %Stublinked(i,:) = outputData;
        %BlinkTimes(i,:) = outputTimes;
        Stublinked = outputData;
        BlinkTimes = outputTimes;
    catch
        %disp(['Pupil data: Stublinks60 crashes: CP ' num2str(i)])
        Stublinked = NaN(1,round(DURATION*SAMPLE_RATE));
        BlinkTimes = ones(1,round(DURATION*SAMPLE_RATE));
        %Stublinked(i,:) = NaN(1,DURATION*SAMPLE_RATE);
        %BlinkTimes(i,:) = ones(1,DURATION*SAMPLE_RATE);
    end
end

end