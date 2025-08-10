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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::inverseDistAndVolDiffusivity

Description
    Calculation of a mixed diffusivity based on inverse volume
    and distance. Two additional coefficent specify:
    powDist - raise the calculated distance to this power
    powVol - raise the calculated volume to this power

SourceFiles
    inverseDistAndVolDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistAndVolDiffusivity_H
#define inverseDistAndVolDiffusivity_H

#include "motionDiffusivity/uniform/uniformDiffusivity.H"
#include "primitives/strings/lists/wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inverseDistAndVolDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class inverseDistAndVolDiffusivity
:
    public uniformDiffusivity
{
    // Private data

        //- Patches selected to base the distance on
        wordReList patchNames_;

        //- Ramp coefficients

            //- raise distance to this power
            scalar powDist_;

            //- raise volume to this power
            scalar powVol_;

    // Private Member Functions

        //- Return patch-cell-centre distance field
        tmp<scalarField> y() const;

        //- Disallow default bitwise copy construct
        inverseDistAndVolDiffusivity(const inverseDistAndVolDiffusivity&);

        //- Disallow default bitwise assignment
        void operator=(const inverseDistAndVolDiffusivity&);


public:

    //- Runtime type information
    TypeName("inverseDistAndVol");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        inverseDistAndVolDiffusivity
        (
            const fvMesh& mesh,
            Istream& mdData
        );


    //- Destructor
    virtual ~inverseDistAndVolDiffusivity();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
