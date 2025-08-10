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
 \*---------------------------------------------------------------------------*/

#include "interpolation/curveFitting/pieceWiseCubicCurveFitter.H"

namespace Foam
{

    pieceWiseCubicCurveFitter::pieceWiseCubicCurveFitter( const vector2DField & data )
    {
        originalData_ = data ;
    }
/*
    void pieceWiseCubicCurveFitter::setInterval( const scalar & input )
    {

        i2_ = -1 ;
        i1_ = -1 ;

        forAll( originalData_ , indexI )
        {
            if (input < originalData_[indexI].x() && i2_ == -1) {
                i2_ = indexI ;
                x2_ = originalData_[i2_].x() ;
                f2_ = originalData_[i2_].y();
                if (i2_ != 0) {
                    i1_ = i2_ - 1 ;
                    x1_ = originalData_[i1_].x() ;
                    f1_ = originalData_[i1_].y() ;
                }
            }
        }
        if (i2_ == -1) {
            i1_ = originalData_.size( ) - 1 ;
        }



    }

    void pieceWiseCubicCurveFitter::setCoefficients( )
    {


        if (i2_ == -1) {
            i2_ = i1_;
            i1_ = i1_ - 1;
            x2_ = originalData_[i2_].x();;
            x1_ = originalData_[i1_].x();
            f2_ = originalData_[i2_].y();
            f1_ = originalData_[i1_].y();
            d1_ = (f2_-f1_)/(x2_-x1_);
            d2_ = d1_;
        }
        else if (i1_ == -1) {
            i1_ = i2_;
            x1_ = originalData_[i1_].x();
            f1_ = f2_;

            i2_ += 1;
            x2_ = originalData_[i2_].x();
            f2_ = originalData_[i2_].y();

            d2_ = (f2_-f1_)/(x2_-x1_);
            d1_ = d2_ ;
        }
        else if (i1_ == 0) {
            d1_ = (-3*originalData_[i1_].y()+4*originalData_[i1_+1].y()-originalData_[i1_+2].y())/2*(originalData_[i1_+1].x()-originalData_[i1_].x());
            d2_ = (originalData_[i2_+1].y()-originalData_[i2_-1].y())/(originalData_[i2_+1].x()-originalData_[i2_-1].x());
        }
        else if (i2_ == originalData_.size()-1){
            d2_ = (3*originalData_[i2_].y()-4*originalData_[i2_-1].y()+originalData_[i2_-2].y())/2*(originalData_[i2_].x()-originalData_[i2_-1].x());
            d1_ = (originalData_[i1_+1].y()-originalData_[i1_-1].y())/(originalData_[i1_+1].x()-originalData_[i1_-1].x());
        }
        else{
            d1_ = (originalData_[i1_+1].y()-originalData_[i1_-1].y())/(originalData_[i1_+1].x()-originalData_[i1_-1].x());
            d2_ = (originalData_[i2_+1].y()-originalData_[i2_-1].y())/(originalData_[i2_+1].x()-originalData_[i2_-1].x());
        }

        A_ = originalData_[i1_].y() ;
        B_ = d1_;
        h_ = x2_-x1_;
        f22_ = f2_-f1_-d1_*h_;
        d22_ = d2_-d1_;
        D_ = (h_*d22_-2*f22_)/pow3(h_);
        C_ = (f22_-D_*pow3(h_))/sqr(h_);

    }
*/
    scalar pieceWiseCubicCurveFitter::interpolate( const scalar & input ) const
    {

        //SET INTERVAL

        scalar i1_,i2_,x1_,x2_,f1_,f2_,d1_,d2_,f22_,d22_,A_,B_,C_,D_,h_;

        i2_ = -1 ;
        i1_ = -1 ;
        x1_ = -1 ;
        x2_ = -1 ;
        f1_ = -1 ;
        f2_ = -1 ;
        d1_ = -1 ;
        d2_ = -1 ;
        f22_ = -1 ;
        d22_ = -1 ;
        A_ = -1 ;
        B_ = -1 ;
        C_ = -1 ;
        D_ = -1 ;
        h_= -1 ;
        forAll( originalData_ , indexI )
        {
            if (input < originalData_[indexI].x() && i2_ == -1) {
                i2_ = indexI ;
                x2_ = originalData_[i2_].x() ;
                f2_ = originalData_[i2_].y();
                if (i2_ != 0) {
                    i1_ = i2_ - 1 ;
                    x1_ = originalData_[i1_].x() ;
                    f1_ = originalData_[i1_].y() ;
                }
            }
        }
        if (i2_ == -1) {
            i1_ = originalData_.size( ) - 1 ;
        }


        //SET COEFFICIENTS

        if (i2_ == -1) {
            i2_ = i1_;
            i1_ = i1_ - 1;
            x2_ = originalData_[i2_].x();;
            x1_ = originalData_[i1_].x();
            f2_ = originalData_[i2_].y();
            f1_ = originalData_[i1_].y();
            d1_ = (f2_-f1_)/(x2_-x1_);
            d2_ = d1_;
        }
        else if (i1_ == -1) {
            i1_ = i2_;
            x1_ = originalData_[i1_].x();
            f1_ = f2_;

            i2_ += 1;
            x2_ = originalData_[i2_].x();
            f2_ = originalData_[i2_].y();

            d2_ = (f2_-f1_)/(x2_-x1_);
            d1_ = d2_ ;
        }
        else if (i1_ == 0) {
            d1_ = (-3*originalData_[i1_].y()+4*originalData_[i1_+1].y()-originalData_[i1_+2].y())/2*(originalData_[i1_+1].x()-originalData_[i1_].x());
            d2_ = (originalData_[i2_+1].y()-originalData_[i2_-1].y())/(originalData_[i2_+1].x()-originalData_[i2_-1].x());
        }
        else if (i2_ == originalData_.size()-1){
            d2_ = (3*originalData_[i2_].y()-4*originalData_[i2_-1].y()+originalData_[i2_-2].y())/2*(originalData_[i2_].x()-originalData_[i2_-1].x());
            d1_ = (originalData_[i1_+1].y()-originalData_[i1_-1].y())/(originalData_[i1_+1].x()-originalData_[i1_-1].x());
        }
        else{
            d1_ = (originalData_[i1_+1].y()-originalData_[i1_-1].y())/(originalData_[i1_+1].x()-originalData_[i1_-1].x());
            d2_ = (originalData_[i2_+1].y()-originalData_[i2_-1].y())/(originalData_[i2_+1].x()-originalData_[i2_-1].x());
        }

        A_ = originalData_[i1_].y() ;
        B_ = d1_;
        h_ = x2_-x1_;
        f22_ = f2_-f1_-d1_*h_;
        d22_ = d2_-d1_;
        D_ = (h_*d22_-2*f22_)/pow3(h_);
        C_ = (f22_-D_*pow3(h_))/sqr(h_);

        return A_ + (input-x1_)*(B_+(input-x1_)*(C_+(input-x1_)*D_));

    }

}
