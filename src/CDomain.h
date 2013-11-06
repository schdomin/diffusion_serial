#ifndef CDOMAIN_H_
#define CDOMAIN_H_

#include <fstream> //ds streaming to file for visualization
#include <math.h>  //ds fabs, etc.
#include <iostream>



namespace Diffusion
{

class CDomain
{

//ds ctor/dtor
public:

    CDomain( const double& p_dDiffusionCoefficient,
             const std::pair< double, double >& p_prBoundaries,
             const double& p_dGridPointSpacing,
             const double& p_dTimeStepSize );
    ~CDomain( );

//ds attributes
private:

    //ds heat structure - we use a standard structure since the constant size is known at allocation and we benefit from access speed
    double** m_gridHeat;

    //ds coefficient matrix
    double** m_matCoefficients;

    //domain properties
    const double m_dDiffusionCoefficient;
    const std::pair< double, double > m_prBoundaries;
    const double m_dGridPointSpacing;
    const unsigned int m_uNumberOfGridPoints1D;
    const double m_dTimeStepSize;
    const double m_dDiffusionFactor;

    //ds stream for offline data - needed for the movie and graphs
    std::string m_strLogHeatDistribution;
    std::string m_strLogNorms;

//ds accessors
public:

    void updateHeatDistributionNumerical( );
    void updateHeatDistributionAnalytical( const double& p_dCurrentTime );
    void saveHeatGridToStream( );
    void saveNormsToStream( const double& p_dCurrentTime );
    void writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const;
    void writeNormsToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize ) const;

//ds helpers
private:

    void setInitialHeatDistribution( );
    double getHeatAnalytical( const double& p_dX, const double& p_dY, const double& p_dT ) const;
    void computeCoefficients( );

    //ds snippet from exercise sheet 1:1
    void solveMatrix( int n, double *a, double *b, double *c, double *v, double *x );

}; //class CDomain

} //namespace Diffusion



#endif //CDOMAIN_H_
