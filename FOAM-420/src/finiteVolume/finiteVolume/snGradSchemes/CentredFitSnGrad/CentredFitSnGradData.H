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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::CentredFitSnGradData

Description
    Data for centred fit snGrad schemes

SourceFiles
    CentredFitSnGradData.C

\*---------------------------------------------------------------------------*/

#ifndef CentredFitSnGradData_H
#define CentredFitSnGradData_H

#include "interpolation/surfaceInterpolation/schemes/FitData/FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedCentredCellToFaceStencil;

/*---------------------------------------------------------------------------*\
                    Class CentredFitSnGradData Declaration
\*---------------------------------------------------------------------------*/

template<class Polynomial>
class CentredFitSnGradData
:
    public FitData
    <
        CentredFitSnGradData<Polynomial>,
        extendedCentredCellToFaceStencil,
        Polynomial
    >
{
    // Private data

        //- For each cell in the mesh store the values which multiply the
        //  values of the stencil to obtain the gradient for each direction
        List<scalarList> coeffs_;


public:

    TypeName("CentredFitSnGradData");


    // Constructors

        //- Construct from components
        CentredFitSnGradData
        (
            const fvMesh& mesh,
            const extendedCentredCellToFaceStencil& stencil,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~CentredFitSnGradData()
    {}


    // Member functions

        //- Return reference to fit coefficients
        const List<scalarList>& coeffs() const
        {
            return coeffs_;
        }

        //- Calculate the fit for the specified face and set the coefficients
        void calcFit
        (
            scalarList& coeffsi, // coefficients to be set
            const List<point>&,  // Stencil points
            const scalar wLin,   // Weight for linear approximation (weights
                                 // nearest neighbours)
            const scalar deltaCoeff, // uncorrected delta coefficient
            const label faci     // Current face index
        );

        void calcFit();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/snGradSchemes/CentredFitSnGrad/CentredFitSnGradData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
