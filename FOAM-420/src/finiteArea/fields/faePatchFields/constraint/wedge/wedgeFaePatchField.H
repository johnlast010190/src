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
    wedgeFaePatchField

Description

SourceFiles
    wedgeFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFaePatchField_H
#define wedgeFaePatchField_H

#include "fields/faePatchFields/faePatchField/faePatchField.H"
#include "faMesh/faPatches/constraint/wedge/wedgeFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class wedgeFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class wedgeFaePatchField
:
    public faePatchField<Type>
{

public:

    //- Runtime type information
    TypeName(wedgeFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        wedgeFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wedgeFaePatchField
        (
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wedgeFaePatchField onto a new patch
        wedgeFaePatchField
        (
            const wedgeFaePatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, faEdgeMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        wedgeFaePatchField
        (
            const wedgeFaePatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return tmp<faePatchField<Type>>
            (
                new wedgeFaePatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wedgeFaePatchField
        (
            const wedgeFaePatchField<Type>&,
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
                new wedgeFaePatchField<Type>(*this, iF)
            );
        }

        //Destructor
        virtual ~wedgeFaePatchField()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "fields/faePatchFields/constraint/wedge/wedgeFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
