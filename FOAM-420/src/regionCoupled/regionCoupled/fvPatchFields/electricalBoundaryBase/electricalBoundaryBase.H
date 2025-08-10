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
    (c) 2019 Esi Ltd.
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::electricalBoundaryBase

Description
    Common functions used in electrical boundaries.

Usage

SourceFiles
    electricalBoundaryBase.C

\*---------------------------------------------------------------------------*/

#ifndef electricalBoundaryBase_H
#define electricalBoundaryBase_H

#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/enums/NamedEnum.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class electricalBoundaryBase Declaration
\*---------------------------------------------------------------------------*/

class electricalBoundaryBase
{
protected:

    // Protected data

        //- Underlying patch
        const fvPatch& patch_;

public:

    // Constructors

        //- Construct from patch
        electricalBoundaryBase
        (
            const fvPatch& patch
        );

        //- Construct from patch and dictionary
        electricalBoundaryBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );

         //- Construct from patch and electricalBoundaryBase
        electricalBoundaryBase
        (
            const fvPatch& patch,
            const electricalBoundaryBase& base
        );


    // Member functions

        //- Return whether anisotropic (tensor conductivity) found or not
        bool anisotropic() const;

        //- Return nf*sigma or nf&sigma in scalar / tensor case respectively
        tmp<vectorField> nfSigma() const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
