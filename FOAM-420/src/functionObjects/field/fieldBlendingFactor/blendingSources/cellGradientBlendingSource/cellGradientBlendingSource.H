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
    (c) 2021 Esi Ltd.

Class
    Foam::cellGradientBlendingSource

Description
    cell centre gradient blending source

SourceFiles
    cellGradientBlendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef cellGradientBlendingSource_H
#define cellGradientBlendingSource_H

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingSources
{

/*---------------------------------------------------------------------------*\
                   Class cellGradientBlendingSource Declaration
\*---------------------------------------------------------------------------*/

class cellGradientBlendingSource
:
    public blendingSource
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellGradientBlendingSource(const cellGradientBlendingSource&);

        //- Disallow default bitwise assignment
        void operator=(const cellGradientBlendingSource&);

public:

    //- Runtime type information
        TypeName("cellGradient");

    // Constructors

        //- Construct from components
        cellGradientBlendingSource(const objectRegistry& obr, const fvMesh& mesh, const dictionary& data);

    //- Destructor
    virtual ~cellGradientBlendingSource(){}


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
