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
                                               m_dTimeStepSize( p_dTimeStepSize )
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
    };

    ~CDomain( )
    {
        //ds deallocate heat structure
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            delete[] m_gridHeat[u];
        }

        delete[] m_gridHeat;
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

    //ds stream for offline data - needed for the movie and graphs
    std::string m_strHeatInformation;

//ds accessors
public:

    void updateHeatDistribution( const double& p_dCurrentTime )
    {
        //ds allocate temporary memory for the data structure (needed for t=1/2)
        double gridHeatNew[m_uNumberOfGridPoints1D][m_uNumberOfGridPoints1D];

        //ds STEP 1: loop over all indexi (i and j are used to conform with exercise parameters)
        for( unsigned int i = 0; i < m_uNumberOfGridPoints1D; ++i )
        {
            //ds structures for implicit equation system (values set to 0 per default)
            double vecNewHeat[m_uNumberOfGridPoints1D];
            double vecPrevHeat[m_uNumberOfGridPoints1D];
            double vecA[m_uNumberOfGridPoints1D];
            double vecB[m_uNumberOfGridPoints1D];
            double vecC[m_uNumberOfGridPoints1D];

            //ds set values of implicit solver vectors
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                vecA[u] =    m_dTimeStepSize*m_dDiffusionCoefficient/( m_dGridPointSpacing*m_dGridPointSpacing );
                vecB[u] =  1-m_dTimeStepSize*m_dDiffusionCoefficient/( m_dGridPointSpacing*m_dGridPointSpacing );
                vecC[u] = -m_dTimeStepSize/2*m_dDiffusionCoefficient/( m_dGridPointSpacing*m_dGridPointSpacing );
            }

            //ds for the current column loop over all rows
            for( unsigned int j = 0; j < m_uNumberOfGridPoints1D; ++j )
            {
                //ds transform coordinates
                const double dX( static_cast< double >( i )/m_uNumberOfGridPoints1D );
                const double dY( static_cast< double >( j )/m_uNumberOfGridPoints1D );

                //ds first add explicit part in y direction
                gridHeatNew[i][j] = m_gridHeat[i][j] + ( m_dDiffusionCoefficient*m_dTimeStepSize/2 )*getHeatDD( dX, dY, p_dCurrentTime );

                //ds and just set previous heat
                vecPrevHeat[j] = m_gridHeat[i][j];
            }

            //ds compute the solution for all current rows
            solveMatrix( m_uNumberOfGridPoints1D, vecA, vecB, vecC, vecPrevHeat, vecNewHeat );

            //ds add the new part in x direction
            for( unsigned int j = 0; j < m_uNumberOfGridPoints1D; ++j )
            {
                //ds directly add to the new heat grid
                gridHeatNew[i][j] += vecNewHeat[j];
            }
        }

        //ds STEP 2: loop over all indexi (i and j are used to conform with exercise parameters)
        for( unsigned int j = 0; j < m_uNumberOfGridPoints1D; ++j )
        {
            //ds structures for implicit equation system (values set to 0 per default)
            double vecNewHeat[m_uNumberOfGridPoints1D];
            double vecPrevHeat[m_uNumberOfGridPoints1D];
            double vecA[m_uNumberOfGridPoints1D];
            double vecB[m_uNumberOfGridPoints1D];
            double vecC[m_uNumberOfGridPoints1D];

            //ds set values of implicit solver vectors
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                vecA[u] =    m_dTimeStepSize*m_dDiffusionCoefficient/( m_dGridPointSpacing*m_dGridPointSpacing );
                vecB[u] =  1-m_dTimeStepSize*m_dDiffusionCoefficient/( m_dGridPointSpacing*m_dGridPointSpacing );
                vecC[u] = -m_dTimeStepSize/2*m_dDiffusionCoefficient/( m_dGridPointSpacing*m_dGridPointSpacing );
            }

            //ds for the current column loop over all rows
            for( unsigned int i = 0; i < m_uNumberOfGridPoints1D; ++i )
            {
                //ds transform coordinates
                const double dX( static_cast< double >( i )/m_uNumberOfGridPoints1D );
                const double dY( static_cast< double >( j )/m_uNumberOfGridPoints1D );

                //ds first add explicit part in x direction
                gridHeatNew[i][j] += ( m_dDiffusionCoefficient*m_dTimeStepSize/2 )*getHeatDD( dX, dY, p_dCurrentTime );

                //ds and just set previous heat
                vecPrevHeat[i] = m_gridHeat[i][j];
            }

            //ds compute the solution for all current columns
            solveMatrix( m_uNumberOfGridPoints1D, vecA, vecB, vecC, vecPrevHeat, vecNewHeat );

            //ds add the new part in y direction
            for( unsigned int i = 0; i < m_uNumberOfGridPoints1D; ++i )
            {
                //ds directly add to the new heat grid
                gridHeatNew[i][j] += vecNewHeat[i];
            }
        }

        //ds STEP 3: update main heat grid
        for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
        {
            for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
            {
                m_gridHeat[u][v] = gridHeatNew[u][v];
            }
        }

        //ds and clear boundaries
        clearBoundaries( );
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
                m_strHeatInformation += chBuffer;
                m_strHeatInformation += " ";
            }

            //ds add new line for new row
            m_strHeatInformation += '\n';
        }
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
            ofsFile << m_uNumberOfGridPoints1D << " " << p_uNumberOfTimeSteps << "\n" << m_strHeatInformation;
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

    //ds doubly derived heat (by x or y)
    double getHeatDD( const double& p_dX, const double& p_dY, const double& p_dT ) const
    {
        //ds formula provided
        return sin( p_dX*2*M_PI )*sin( p_dY*2*M_PI )*exp( -8*m_dDiffusionCoefficient*M_PI*M_PI*p_dT )*( -4*M_PI*M_PI );
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
