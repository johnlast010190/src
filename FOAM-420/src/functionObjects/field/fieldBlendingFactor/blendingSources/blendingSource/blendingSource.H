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
    Foam::blendingSource

Description
    Abstract base class for blending functions

SourceFiles
    blendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef blendingSource_H
#define blendingSource_H

#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                          Class blendingSource Declaration
\*---------------------------------------------------------------------------*/

class blendingSource
{
protected:

    // Protected Data

        //- value to coefficient translator
        autoPtr<blendingFunction> ramp_;

private:

    // Private Data

        //- output switch
        Switch writeBlendingField_;

    // Protected Member Functions

        //- calculate and apply blending factor
        tmp<surfaceScalarField> convertValueToBlendingFactor
        (
            const surfaceScalarField& valueField,
            bool verbose = true
        );

protected:

    // Protected Data
        //- object registry
        const objectRegistry& obr_;

        //- const reference to mesh
        const fvMesh& mesh_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        blendingSource(const blendingSource&);

        //- Disallow default bitwise assignment
        void operator=(const blendingSource&);


public:

    //- Runtime type information
        TypeName("blendingSource");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            blendingSource,
            dictionary,
            (
                const objectRegistry& obr,
                const fvMesh& mesh,
                const dictionary& blendingSourceDict
            ),
            (obr, mesh, blendingSourceDict)
        );

    // Constructors

        //- Construct from components
        blendingSource(const objectRegistry& obr, const fvMesh& mesh, const dictionary& data);

    // Selectors

        //- Return an autoPtr to the selected blendingSource
        static autoPtr<blendingSource> New
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& blendingSourceDict
        );


    //- Destructor
    virtual ~blendingSource()
    {}


    // Member Functions

        //- return blending field from this
        tmp<surfaceScalarField> blendingFactor();

        // write blending source field to disk
        void write();

        //- calculate source field
        virtual tmp<surfaceScalarField> sourceField() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
