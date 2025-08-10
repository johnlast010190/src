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
    (c) 2017 Esi Ltd.

Class
    Foam::phaseBlendingSource

Description
    Blending source for multiphase interface. Indicator field, 1 at interface,
    0 elsewhere

SourceFiles
    phaseBlendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef phaseBlendingSource_H
#define phaseBlendingSource_H

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingSources
{

/*---------------------------------------------------------------------------*\
                   Class phaseBlendingSource Declaration
\*---------------------------------------------------------------------------*/

class phaseBlendingSource
:
    public blendingSource
{
private:

    // Private Data



    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseBlendingSource(const phaseBlendingSource&);

        //- Disallow default bitwise assignment
        void operator=(const phaseBlendingSource&);

public:

    //- Runtime type information
        TypeName("phase");

    // Constructors

        //- Construct from components
        phaseBlendingSource(const objectRegistry& obr, const fvMesh& mesh, const dictionary& data);

    //- Destructor
    virtual ~phaseBlendingSource(){}


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
