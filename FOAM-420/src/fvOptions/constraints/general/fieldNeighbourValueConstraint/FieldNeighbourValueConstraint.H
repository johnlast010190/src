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
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::fv::FieldNeighbourValueConstraint

Group
    grpFvOptionsConstraints

Description
    Constrain the field values within a specified region based on current average neighbour values.

Usage
    For example to remove bad cells from matrix solve:
    \verbatim
    cellRemoval
    {
        type            scalarFieldNeighbourValueConstraint;
        active          yes;

        selectionMode   cellZone;
        cellZone        cellRemoval;
        fieldNames
        {
            k         ;
            epsilon   ;
            T         ;
            p_rgh     ;
        }
    }
    \endverbatim

See also
    Foam::fvOption

SourceFiles
    FieldNeighbourValueConstraint.C
    FieldNeighbourValueConstraints.C

\*---------------------------------------------------------------------------*/

#ifndef FieldNeighbourValueConstraint_H
#define FieldNeighbourValueConstraint_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class FieldNeighbourValueConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FieldNeighbourValueConstraint
:
    public cellSetOption
{
public:

    //- Runtime type information
    TypeName("FieldNeighbourValueConstraint");


    // Constructors

        //- Construct from components
        FieldNeighbourValueConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Set value on field
        virtual void constrain(fvMatrix<Type>& eqn, const label fieldi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FieldNeighbourValueConstraint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
