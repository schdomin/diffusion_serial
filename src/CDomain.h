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
             const double& p_dTimeStepSize ) : m_dDiffusionCoefficient( p_dDiffusionCoefficient ),
                                               m_prBoundaries( p_prBoundaries ),
                                               m_dGridPointSpacing( p_dGridPointSpacing ),
                                               m_uNumberOfGridPoints1D( ceil( fabs( p_prBoundaries.first )+fabs( p_prBoundaries.second ) )/m_dGridPointSpacing+1 ),
                                               m_dTimeStepSize( p_dTimeStepSize ),
                                               m_dDiffusionFactor( p_dTimeStepSize*p_dDiffusionCoefficient/( 2*( p_dGridPointSpacing*p_dGridPointSpacing ) ) ),
                                               m_strLogHeatDistribution( "" ),
                                               m_strLogTotalHeat( "" )

    {
        //ds allocate memory for the data structure
        m_gridHeat = new double*[m_uNumberOfGridPoints1D];

        //ds for each element
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            m_gridHeat[u] = new double[m_uNumberOfGridPoints1D];
        }

        //ds initialize grid
        setInitialHeatDistribution( );

        //ds allocate memory for the support structure
        m_matCoefficients = new double*[3];

        //ds and for the 3 diagonals
        m_matCoefficients[0] = new double[m_uNumberOfGridPoints1D-1]; //sup
        m_matCoefficients[1] = new double[m_uNumberOfGridPoints1D];   //main
        m_matCoefficients[2] = new double[m_uNumberOfGridPoints1D-1]; //sub
    };

    ~CDomain( )
    {
        //ds deallocate heat structure
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            delete[] m_gridHeat[u];
        }

        delete[] m_gridHeat;

        //ds deallocate coefficient structure
        delete[] m_matCoefficients[0];
        delete[] m_matCoefficients[1];
        delete[] m_matCoefficients[2];
        delete[] m_matCoefficients;
    };

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
    std::string m_strLogTotalHeat;

//ds accessors
public:

    void updateHeatDistributionNumerical( )
    {
        //ds STEP 1: loop over all rows (i and j are used to conform with exercise parameters)
        for( unsigned int j = 0; j < m_uNumberOfGridPoints1D; ++j )
        {
            //ds structures for implicit equation system (values set to 0 per default)
            double vecNewHeat[m_uNumberOfGridPoints1D];
            double vecPreviousHeat[m_uNumberOfGridPoints1D];

            //ds boundaries
            vecPreviousHeat[0]                         = 0.0;
            vecPreviousHeat[m_uNumberOfGridPoints1D-1] = 0.0;

            //ds for the current row get the current RHS
            for( unsigned int i = 1; i < m_uNumberOfGridPoints1D-1; ++i )
            {
                //ds get previous heat
                vecPreviousHeat[i] = m_gridHeat[i][j] + m_dDiffusionFactor*( m_gridHeat[i+1][j] - 2*m_gridHeat[i][j] + m_gridHeat[i-1][j] );
            }

            //ds update coefficient matrix
            computeCoefficients( );

            //ds compute the new heat
            solveMatrix( m_uNumberOfGridPoints1D, m_matCoefficients[0], m_matCoefficients[1], m_matCoefficients[2], vecPreviousHeat, vecNewHeat );

            //ds set the values to the grid
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                m_gridHeat[u][j] = vecNewHeat[u];
            }
        }

        //ds STEP 2: loop over all columns (i and j are used to conform with exercise parameters)
        for( unsigned int i = 0; i < m_uNumberOfGridPoints1D; ++i )
        {
            //ds structures for implicit equation system (values set to 0 per default)
            double vecNewHeat[m_uNumberOfGridPoints1D];
            double vecPreviousHeat[m_uNumberOfGridPoints1D];

            //ds boundaries
            vecPreviousHeat[0]                         = 0.0;
            vecPreviousHeat[m_uNumberOfGridPoints1D-1] = 0.0;

            //ds for the current row get the current RHS
            for( unsigned int j = 1; j < m_uNumberOfGridPoints1D-1; ++j )
            {
                //ds get previous heat
                vecPreviousHeat[j] = m_gridHeat[i][j] + m_dDiffusionFactor*( m_gridHeat[i][j+1] - 2*m_gridHeat[i][j] + m_gridHeat[i][j-1] );
            }

            //ds updates coefficient matrix
            computeCoefficients( );

            //ds compute the new heat
            solveMatrix( m_uNumberOfGridPoints1D, m_matCoefficients[0], m_matCoefficients[1], m_matCoefficients[2], vecPreviousHeat, vecNewHeat );

            //ds set the values to the grid
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                m_gridHeat[i][u] = vecNewHeat[u];
            }
        }
    }

    void updateHeatDistributionAnalytical( const double& p_dCurrentTime )
    {
        /*ds for all grid points
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
            {
                //ds get the exact heat at this point
                m_gridHeat[u][v] = sin( u*m_dGridPointSpacing*2*M_PI )*sin( v*m_dGridPointSpacing*2*M_PI )*exp( -8*m_dDiffusionCoefficient*M_PI*M_PI*p_dCurrentTime );
            }
        }*/
    }

    void saveHeatGridToStream( )
    {
        //ds buffer for snprintf
        char chBuffer[16];

        //ds add each element
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
            {
                //ds get the integrals stream
                std::snprintf( chBuffer, 16, "%f", m_gridHeat[u][v] );

                //ds add buffer and space to string
                m_strLogHeatDistribution += chBuffer;
                m_strLogHeatDistribution += " ";
            }

            //ds add new line for new row
            m_strLogHeatDistribution += '\n';
        }
    }

    void saveTotalHeatToStream( )
    {
        //ds get total heat
        double dTotalHeat( 0.0 );

        //ds for all grid points
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
            {
                //ds add the current heat value
                dTotalHeat += m_gridHeat[u][v];
            }
        }

        //ds buffer for snprintf
        char chBuffer[16];

        //ds get the integrals stream
        std::snprintf( chBuffer, 16, "%f", dTotalHeat );

        //ds append the buffer to our string
        m_strLogTotalHeat += chBuffer;
        m_strLogTotalHeat += "\n";
    }

    void writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const
    {
        //ds ofstream object
        std::ofstream ofsFile;

        //ds open the file for writing
        ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

        //ds if it worked
        if( ofsFile.is_open( ) )
        {
            //ds first dump setup information number of points and time steps
            ofsFile << m_uNumberOfGridPoints1D << " " << p_uNumberOfTimeSteps << "\n" << m_strLogHeatDistribution;
        }

        //ds close the file
        ofsFile.close( );
    }

    void writeTotalHeatToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize )
    {
        //ds ofstream object
        std::ofstream ofsFile;

        //ds open the file for writing
        ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

        //ds if it worked
        if( ofsFile.is_open( ) )
        {
            //ds dump information to file
            ofsFile << p_uNumberOfTimeSteps << " " << p_dTimeStepSize << "\n" << m_strLogTotalHeat;
        }

        //ds close the file
        ofsFile.close( );
    }

//ds helpers
private:

    void setInitialHeatDistribution( )
    {
        //ds loop over all indexi
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
            {
                //ds set initial heat value
                m_gridHeat[u][v] = sin( u*m_dGridPointSpacing*2*M_PI )*sin( v*m_dGridPointSpacing*2*M_PI );
            }
        }

        //ds set boundary values
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            m_gridHeat[u][0]                         = 0.0;
            m_gridHeat[u][m_uNumberOfGridPoints1D-1] = 0.0;
            m_gridHeat[0][u]                         = 0.0;
            m_gridHeat[m_uNumberOfGridPoints1D-1][u] = 0.0;
        }
    }

    void computeCoefficients( )
    {
        //ds set the coefficients
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            //ds sup, main and sub diagonal
            m_matCoefficients[0][u] = -m_dDiffusionFactor;
            m_matCoefficients[1][u] = 1+2*m_dDiffusionFactor;
            m_matCoefficients[2][u] = -m_dDiffusionFactor;
        }

        //ds sup diagonal: 1 .. n-1 size: n-1
        m_matCoefficients[0][0]                         = 0.0;
        m_matCoefficients[0][m_uNumberOfGridPoints1D-1] = 0.0;

        //main diagonal boundary condition
        m_matCoefficients[1][0]                         = 1.0;
        m_matCoefficients[1][m_uNumberOfGridPoints1D-1] = 1.0;

        //ds sub diagonal: 0 .. n-2 size: n-1
        m_matCoefficients[2][0]                         = 0.0;
        m_matCoefficients[2][m_uNumberOfGridPoints1D-1] = 0.0;
    }

    //ds snippet from exercise sheet 1:1
    void solveMatrix( int n, double *a, double *b, double *c, double *v, double *x )
    {
        /*
        * n − number of equations
        * a − sub−diagonal (means it is the diagonal below the main diagonal) −− indexed from 1..n−1
        * b − the main diagonal
        * c − sup−diagonal (means it is the diagonal above the main diagonal) −− indexed from 0..n−2
        * v − right part
        * x − the answer
        */

        for( int i = 1; i < n; i++ )
        {
            double m = a[i] / b[i-1];
            b[i]     = b[i] - m*c[i-1];
            v[i]     = v[i] - m*v[i-1];
        }

        x[n-1] = v[n-1] / b[n-1];

        for( int i = n-2; i >= 0; i-- )
        {
            x[i] = ( v[i]-c[i] * x[i+1] )/b[i];
        }
    }

}; //class CDomain

} //namespace Diffusion



#endif //CDOMAIN_H_
