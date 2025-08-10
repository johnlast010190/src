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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::nonConformalProcessorCyclicPointPatchField

Description
    Constraint pointPatchField for nonConformalProcessorCyclic couplings.

SourceFiles
    nonConformalProcessorCyclicPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicPointPatchField_H
#define nonConformalProcessorCyclicPointPatchField_H

#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "nonConformal/pointPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nonConformalProcessorCyclicPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicPointPatchField
:
    public pointPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPointPatch::typeName_());


    // Constructors

        //- Construct from patch and nonConformalProcessorCyclic field
        nonConformalProcessorCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, nonConformalProcessorCyclic field and
        //  dictionary
        nonConformalProcessorCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalProcessorCyclicPointPatchField
        (
            const nonConformalProcessorCyclicPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting nonConformalProcessorCyclic field
        //  reference
        nonConformalProcessorCyclicPointPatchField
        (
            const nonConformalProcessorCyclicPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalProcessorCyclicPointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct and return a clone setting nonConformalProcessorCyclic
        //  field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const override
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalProcessorCyclicPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformal/pointPatchFields/nonConformalProcessorCyclic/nonConformalProcessorCyclicPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
