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
    Foam::fieldBlendingItem

Description
    Management class for variable specific blending field

    @verbatim
        {
            stabilise 10;

            blendingSources
            (
                cfl
                {
                    minValue 2;
                    maxValue 5;
                    minBlend 1;
                    maxBlend 0;
                    rampType cubic;
                    write   false;
                }

                Pe
                {
                    minValue 10;
                    maxValue 20;
                    minBlend 1;
                    maxBlend 0.9;
                    rampType cubic;
                }
            );
        }
    @endverbatim

SourceFiles
    fieldBlendingItem.C

\*---------------------------------------------------------------------------*/

#ifndef fieldBlendingItem_H
#define fieldBlendingItem_H

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendingSource;

/*---------------------------------------------------------------------------*\
                      Class fieldBlendingItem Declaration
\*---------------------------------------------------------------------------*/

class fieldBlendingItem
{
private:

    // Private Data
        //- objecti registry
        const objectRegistry& obr_;

        //- const reference to mesh
        const fvMesh& mesh_;

        //- field name
        word fieldName_;

        //- stabilisation coefficient
        scalar stab_;

        //- dynamic blending sources
        PtrList<blendingSource> blendSources_;

    // Private Member Functions

        //- copy, stabilise and limit new blending factor
        void blendingFactorTransfer(scalar& oldBF, scalar newBF);

        //- Disallow default bitwise copy construct
        fieldBlendingItem(const fieldBlendingItem&);

        //- Disallow default bitwise assignment
        void operator=(const fieldBlendingItem&);

public:

    //- Runtime type information
    TypeName("fieldBlendingItem");

    // Constructors

        //- Construct from mesh, base field name and input dictionary
        fieldBlendingItem
        (
            const objectRegistry& obr,
            const fvMesh& mesh,
            const dictionary& dict,
            const word blendingFieldName
        );

    // Destructor

        virtual ~fieldBlendingItem();


    // Member Functions

        // Access

            //- return field name
            word name() const
            {
                return fieldName_;
            }

            //- return const access to named field
            inline const surfaceScalarField& field() const
            {
                return obr_.lookupObject<surfaceScalarField>(fieldName_);
            }


        // Edit

            //- update fieldBlendingItem if dictionary has changed
            void readBlendingItem
            (
                const dictionary& d
            );

        // Update

            //- update field blending factors
            void update(const surfaceScalarField* bsbfPtr = nullptr);

            //- write field blending factors
            void write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
