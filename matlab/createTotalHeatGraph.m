function createIntegralsGraphs( )

%ds filepath
strFilepath = '../bin/totalheat.txt';

%ds open the file
fileID = fopen( strFilepath );

%ds get the first line
cCell = textscan( fileID, '%u %f', 1 );

%ds get timesteps T
uNumberOfTimesteps = cCell{1};
dTimestepSize      = cCell{2};

%ds informative
disp( ['Number of timesteps: ', num2str( uNumberOfTimesteps ) ] );
disp( ['Timestep size: ', num2str( dTimestepSize ) ] );
disp( [ 'starting data import from: ', strFilepath ] ); 
tic;

%ds get the remaining lines
cCell = textscan( fileID, '%f', uNumberOfTimesteps );

vecTotalHeat = cCell{:};

disp( [ 'finished data import - time: ', num2str( toc ) ] );

%ds allocate the timeline
vecTimeline = zeros( uNumberOfTimesteps, 1 );

%ds current time for the vector (safety for integer/double multiplicatioo)
dCurrentTime = 0.0;

%ds fill the timeline
for u = 1:1:uNumberOfTimesteps
    
    %ds set the element
    vecTimeline( u ) = dCurrentTime;
    
    %ds update current time
    dCurrentTime = dCurrentTime+dTimestepSize;
    
end


hFigure1 = figure( 1 );
plot( vecTimeline, vecTotalHeat );
%axis([0, dMaxTime, 0, 1.5*dMaxEnergy]);
title( 'Total Heat' );
xlabel( 'Time' );
ylabel( 'Heat' );
legend( 'Total Heat' );

disp( 'exporting figures as jpg' );

saveas( hFigure1, 'totalheat.jpg' );

disp( 'function ended successfully' );

end
