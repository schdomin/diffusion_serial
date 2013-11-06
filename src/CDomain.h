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
        setHeatDistribution( );

        //ds allocate memory for the support structure
        m_matCoefficients = new double*[3];

        //ds for the 3 diagonals
        m_matCoefficients[0] = new double[m_uNumberOfGridPoints1D-1];
        m_matCoefficients[1] = new double[m_uNumberOfGridPoints1D];
        m_matCoefficients[2] = new double[m_uNumberOfGridPoints1D-1];

        computeCoefficients( );
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

    //domain properties
    const double m_dDiffusionCoefficient;
    const std::pair< double, double > m_prBoundaries;
    const double m_dGridPointSpacing;
    const unsigned int m_uNumberOfGridPoints1D;
    const double m_dTimeStepSize;

    //ds coefficient matrix
    double** m_matCoefficients;
    const double m_dDiffusionFactor;

    //ds stream for offline data - needed for the movie and graphs
    std::string m_strLogHeatDistribution;
    std::string m_strLogTotalHeat;

//ds accessors
public:

    void updateHeatDistribution( )
    {
        //ds STEP 1: loop over all rows (i and j are used to conform with exercise parameters)
        for( unsigned int j = 0; j < m_uNumberOfGridPoints1D; ++j )
        {
            //ds structures for implicit equation system (values set to 0 per default)
            double vecNewHeat[m_uNumberOfGridPoints1D];
            double vecPreviousHeat[m_uNumberOfGridPoints1D];

            //ds updates coefficient matrix
            computeCoefficients( );

            //ds boundaries
            vecPreviousHeat[0]                         = 0.0;
            vecPreviousHeat[m_uNumberOfGridPoints1D-1] = 0.0;

            //ds for the current row get the current RHS
            for( unsigned int i = 1; i < m_uNumberOfGridPoints1D-1; ++i )
            {
                //ds get previous heat
                vecPreviousHeat[i] = m_gridHeat[i][j] + m_dDiffusionFactor*( m_gridHeat[i-1][j] - 2*m_gridHeat[i][j] + m_gridHeat[i+1][j] );
            }

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

            //ds updates coefficient matrix
            computeCoefficients( );

            //ds boundaries
            vecPreviousHeat[0]                         = 0.0;
            vecPreviousHeat[m_uNumberOfGridPoints1D-1] = 0.0;

            //ds for the current row get the current RHS
            for( unsigned int j = 1; j < m_uNumberOfGridPoints1D-1; ++j )
            {
                //ds get previous heat
                vecPreviousHeat[j] = m_gridHeat[i][j] + m_dDiffusionFactor*( m_gridHeat[i][j-1] - 2*m_gridHeat[i][j] + m_gridHeat[i][j+1] );
            }

            //ds compute the new heat
            solveMatrix( m_uNumberOfGridPoints1D, m_matCoefficients[0], m_matCoefficients[1], m_matCoefficients[2], vecPreviousHeat, vecNewHeat );

            //ds set the values to the grid
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                m_gridHeat[i][u] = vecNewHeat[u];
            }
        }
    }

    void saveHeatGridToStream( )
    {
        //ds buffer for snprintf
        char chBuffer[16];

        //ds add each element to the string
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

    void setHeatDistribution( const double& p_dT = 0.0 )
    {
        //ds loop over all indexi
        for( unsigned int uX = 0; uX < m_uNumberOfGridPoints1D; ++uX )
        {
            for( unsigned int uY = 0; uY < m_uNumberOfGridPoints1D; ++uY )
            {
                //ds transform coordinates
                const double dX( static_cast< double >( uX )/m_uNumberOfGridPoints1D );
                const double dY( static_cast< double >( uY )/m_uNumberOfGridPoints1D );

                //ds set heat value
                m_gridHeat[uX][uY] = getHeat( dX, dY, p_dT );
            }
        }

        //ds and clean up borders
        clearBoundaries( );
    }

    double getHeat( const double& p_dX, const double& p_dY, const double& p_dT ) const
    {
        //ds formula provided
        return sin( p_dX*2*M_PI )*sin( p_dY*2*M_PI )*exp( -8*m_dDiffusionCoefficient*M_PI*M_PI*p_dT );
    }

    //ds boundary conditions
    void clearBoundaries( )
    {
        //ds clear all trivial rows
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
        for( unsigned int u = 1; u < m_uNumberOfGridPoints1D-1; ++u )
        {
            //ds lower, main and upper diagonal
            m_matCoefficients[0][u] = -m_dDiffusionFactor;
            m_matCoefficients[1][u] = 1+2*m_dDiffusionFactor;
            m_matCoefficients[2][u] = -m_dDiffusionFactor;
        }

        //ds additional element for main diagonal
        m_matCoefficients[1][m_uNumberOfGridPoints1D-2] = 1+2*m_dDiffusionFactor;

        //ds boundaries
        m_matCoefficients[0][0]                         = 0.0;
        m_matCoefficients[1][0]                         = 1.0;
        m_matCoefficients[2][0]                         = 0.0;
        m_matCoefficients[0][m_uNumberOfGridPoints1D-2] = 0.0;
        m_matCoefficients[1][m_uNumberOfGridPoints1D-1] = 1.0;
        m_matCoefficients[2][m_uNumberOfGridPoints1D-2] = 0.0;
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
