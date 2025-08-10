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

Class
    Foam::fv::FixedValueConstraint

Group
    grpFvOptionsConstraints

Description
    Constrain the field values within a specified region.

Usage
    For example to set the turbulence properties within a porous region:
    \verbatim
    porosityTurbulence
    {
        type            scalarFixedValueConstraint;
        active          yes;

        selectionMode   cellZone;
        cellZone        porosity;
        fieldValues
        {
            k           1;
            epsilon     150;
        }
    }
    \endverbatim

See also
    Foam::fvOption

SourceFiles
    FixedValueConstraint.C
    fixedValueConstraints.C

\*---------------------------------------------------------------------------*/

#ifndef FixedValueConstraint_H
#define FixedValueConstraint_H

#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class FixedValueConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FixedValueConstraint
:
    public cellSetOption
{
    // Private member data

        //- Field values
        List<Type> fieldValues_;


public:

    //- Runtime type information
    TypeName("FixedValueConstraint");


    // Constructors

        //- Construct from components
        FixedValueConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary&);

        //- Set value on field
        virtual void constrain(fvMatrix<Type>&, const label);

        //- Set value on field
        virtual void constrain(fvBlockMatrix<Type>&, const label);
};

// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

//- Set value on field
template<>
void FixedValueConstraint<vector>::constrain
(
    fvBlockMatrix<vector>&,
    const label
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FixedValueConstraint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
