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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Class
    processorFaePatchField

Description

SourceFiles
    processorFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef processorFaePatchField_H
#define processorFaePatchField_H

#include "fields/faePatchFields/basic/coupled/coupledFaePatchField.H"
#include "faMesh/faPatches/constraint/processor/processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class processorFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorFaePatchField
:
    public coupledFaePatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const processorFaPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(processorFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        processorFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct from patch and internal field and patch field
        processorFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        processorFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given processorFaePatchField onto a new patch
        processorFaePatchField
        (
            const processorFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        processorFaePatchField(const processorFaePatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return tmp<faePatchField<Type>>
            (
                new processorFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        processorFaePatchField
        (
            const processorFaePatchField<Type>&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, faEdgeMesh>& iF
        ) const
        {
            return tmp<faePatchField<Type>>
            (
                new processorFaePatchField<Type>(*this, iF)
            );
        }


    // Destructor

        virtual ~processorFaePatchField();


    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const
            {
                if (Pstream::parRun())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fields/faePatchFields/constraint/processor/processorFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
