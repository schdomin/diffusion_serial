#include "CDomain.h" //ds grid
#include "Timer.h"   //ds time measurement

int main( int argc, char** argv )
{
    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds get parameters
    const double dDiffusionCoefficient( 1.0 );
    const std::pair< double, double > prBoundaries( 0.0, 1.0 );
    const double dGridPointSpacing( 0.01 );
    const unsigned int uNumberOfTimeSteps( 1000 );
    const double dTimeStepSize( 0.5*dGridPointSpacing*dGridPointSpacing/dDiffusionCoefficient );

    //ds user information
    std::cout << "\n---------------------------- DIFFUSION SERIAL SETUP ----------------------------" << std::endl;
    std::cout << "Diffusion Coefficient: " << dDiffusionCoefficient << std::endl;
    std::cout << "        Boundary (2D): [" << prBoundaries.first << ", " << prBoundaries.second << "]" << std::endl;
    std::cout << "   Grid Point Spacing: " << dGridPointSpacing << std::endl;
    std::cout << " Number of Time Steps: " << uNumberOfTimeSteps << std::endl;
    std::cout << "       Time Step Size: " << dTimeStepSize << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;

    //ds allocate domain (automatically creates initial density distribution)
    Diffusion::CDomain cDomain( dDiffusionCoefficient, prBoundaries, dGridPointSpacing, dTimeStepSize );

    //ds information
    std::cout << "               Status:  0% done - current step: 0";

    //ds start simulation
    for( unsigned int uCurrentTimeStep = 1; uCurrentTimeStep < uNumberOfTimeSteps+1; ++uCurrentTimeStep )
    {
        //ds calculate percentage done
        const double dPercentageDone( 100.0*uCurrentTimeStep/uNumberOfTimeSteps );

        //ds get a formatted string -> 100% -> 3 digits
        char chBuffer[4];

        //ds fill the buffer
        std::snprintf( chBuffer, 4, "%3.0f", dPercentageDone );

        //ds print info
        std::cout << '\xd';
        std::cout << "               Status: " << chBuffer << "% done - current step: " << uCurrentTimeStep;

        //ds current time
        const double dCurrentTime( uCurrentTimeStep*dTimeStepSize );

        //ds update domain
        cDomain.updateHeatDistribution( dCurrentTime );

        //ds streaming
        cDomain.saveHeatGridToStream( );
        cDomain.saveTotalHeatToStream( );
    }

    //ds save the stream to a file
    cDomain.writeHeatGridToFile( "bin/simulation.txt", uNumberOfTimeSteps );
    cDomain.writeTotalHeatToFile( "bin/totalheat.txt", uNumberOfTimeSteps, dTimeStepSize );

    //ds stop timing
    const double dDurationSeconds( tmTimer.stop( ) );

    //ds cause an output ostream
    std::cout << std::endl;
    std::cout << "     Computation time: " << dDurationSeconds << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    return 0;
}
