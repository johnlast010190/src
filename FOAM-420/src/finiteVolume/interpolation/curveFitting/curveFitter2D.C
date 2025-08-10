/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2010-2012 Esi Ltd.

 @author Jacques Papper
 \*---------------------------------------------------------------------------*/
#include <vector>

#include "interpolation/curveFitting/curveFitter2D.H"
#include "interpolation/curveFitting/pieceWiseCubicCurveFitter.H"


namespace Foam
{

    curveFitter2D::curveFitter2D( )
    {
    }

    curveFitter2D::curveFitter2D( const List< Tuple2<scalar , vector2DField >> & allCurves )
    {

        forAll( allCurves , indexI )
        {
            pieceWiseCubicCurveFitter curve( allCurves[indexI].second() ) ;
            allDataX_.push_back( allCurves[indexI].first()) ;
            alldata_.push_back( curve ) ;
        }
    }

    curveFitter2D::curveFitter2D( const curveFitter2D& orig )
    {
        this->allDataX_ = orig.allDataX_ ;
        this->alldata_ = orig.alldata_ ;
    }

    curveFitter2D::~curveFitter2D( )
    {
    }

    scalar curveFitter2D::interpolate( const scalar& input1, const scalar& input2 )
    {
        vector2DField newCurve( allDataX_.size( ) ) ;

        forAll( newCurve , indexI )
        {
            newCurve[indexI].x( ) = allDataX_[indexI] ;
            newCurve[indexI].y( ) = alldata_[indexI].interpolate( input1 ) ;
        }
        pieceWiseCubicCurveFitter newInterpolator( newCurve ) ;
        return newInterpolator.interpolate( input2 ) ;
    }

    void curveFitter2D::addAllCurves(const List< Tuple2<scalar,vector2DField>> & allCurves)
    {
        forAll( allCurves , indexI )
        {
            pieceWiseCubicCurveFitter curve( allCurves[indexI].second() ) ;
            allDataX_.push_back( allCurves[indexI].first()) ;
            alldata_.push_back( curve ) ;
        }
    }

    void curveFitter2D::addCurve( const scalar & X , const vector2DField & curveData )
    {

        pieceWiseCubicCurveFitter curve( curveData ) ;
        allDataX_.push_back( X ) ;
        alldata_.push_back( curve ) ;

    }

}
