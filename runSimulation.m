function runSimulation(beta, gridSize, dt, endTime)
    U = zeros(gridSize, gridSize);
    V = zeros(gridSize, gridSize);
    U(floor(gridSize/2), floor(gridSize/2)) = 1;

    video = VideoWriter(sprintf('Euler_beta_%g.mp4', beta), 'MPEG-4');
    video.FrameRate = 20; % frame rate; can be adjusted
    open(video);
    
    frames = cell(1, round(endTime/dt));
    frameIndex = 1;
    
    for t = 0:dt:endTime
        % [U, V] = latticeStepRK4(U, V, dt, beta);
        [U, V] = latticeStep(U, V, dt, beta);
        frames{frameIndex} = U;
        frameIndex = frameIndex + 1;
    end

    for i = 1:length(frames)
        imagesc(frames{i});
        colorbar;
        caxis([-1, 1]);
        title(sprintf('Wave Propagation at Beta = %g, dt = %g, Time = %g', beta, dt, i*dt));
        
        % for video
        frame = getframe(gcf);
        writeVideo(video, frame);
        
        pause(0.0005); % speed of vid playing 
    end

    close(video);  % Close the video file
end
