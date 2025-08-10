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
    (c) 2011 Esi Ltd.

Class
    Foam::zonalDESBlendingSource

Description
    zonalDES blending source

SourceFiles
    zonalDESBlendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef zonalDESBlendingSource_H
#define zonalDESBlendingSource_H

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingSources
{

/*---------------------------------------------------------------------------*\
                   Class zonalDESBlendingSource Declaration
\*---------------------------------------------------------------------------*/

class zonalDESBlendingSource
:
    public blendingSource
{
private:

    // Private Data

        //- CDES length scale constant
        scalar CDES_;

        //- smoothing iterations
        label nSmooth_;

        //- zone indicator field
        autoPtr<surfaceScalarField> fZone_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        zonalDESBlendingSource(const zonalDESBlendingSource&);

        //- Disallow default bitwise assignment
        void operator=(const zonalDESBlendingSource&);

public:

    //- Runtime type information
        TypeName("zonalDES");

    // Constructors

        //- Construct from components
        zonalDESBlendingSource(const objectRegistry& obr, const fvMesh& mesh, const dictionary& data);

    //- Destructor
    virtual ~zonalDESBlendingSource(){}


    // Member Functions

        //- source field calculation
        virtual tmp<surfaceScalarField> sourceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingSources
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
