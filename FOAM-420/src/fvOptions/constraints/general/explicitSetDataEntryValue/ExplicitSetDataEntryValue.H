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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::fv::explicitSetDataEntryValue

Description
    Set values field values explicity.

    Sources described by:

        \verbatim
        <Type>ExplicitSetDataEntryValueCoeffs
        {
            sourceData
            {
                p table
                (
                    (0 0)
                    (1 10)
                    (2 0)
                );
                epsilon  constant   1.5;
            }
        }
        \endverbatim

SeeAlso
    Foam::fvOption

SourceFiles
    explicitSetDataEntryValue.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitSetDataEntryValue_H
#define ExplicitSetDataEntryValue_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "primitives/Tuple2/Tuple2.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class explicitSetDataEntryValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ExplicitSetDataEntryValue
:
    public cellSetOption
{

protected:

    // Protected data

        //- Source value per field
        PtrList<Function1<Type>> sourceData_;


    // Protected functions

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("explicitSetDataEntryValue");


    // Constructors

        //- Construct from components
        ExplicitSetDataEntryValue
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        // Evaluation

            //- Set value on field
            virtual void setValue(fvMatrix<Type>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ExplicitSetDataEntryValue.C"
#   include "ExplicitSetDataEntryValueIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
