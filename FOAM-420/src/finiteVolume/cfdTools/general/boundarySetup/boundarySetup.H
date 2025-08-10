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
    (c) 2010-2016 Esi Ltd.

Class
    Foam::boundarySetup

Description
    Interpret a structured input dictionary to create a boundary for
    a volField. Note this is not the standard boundaryField format dictionary.

SourceFiles
    boundarySetup.C

\*---------------------------------------------------------------------------*/

#ifndef boundarySetup_H
#define boundarySetup_H

#include "fvMesh/fvMesh.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boundarySetup Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class boundarySetup
:
    public PtrList<fvPatchField<Type>>
{
    // boundaryField dictionary for CBC usage
    autoPtr<dictionary> dictPtr_;

    //- default fvPatch-type wide boundaries
    void defaultBoundaries
        (const fvMesh&, const DimensionedField<Type, volMesh>&);

    //- default fvPatch-type wide boundary dictionary entries
    void defaultBoundaries
        (const fvMesh&);

    //- fvPatch-type wide boundaries
    void typeBoundaries
    (
        const fvMesh&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );

    //- fvPatch-type wide boundary dictionary entries
    void typeBoundaries
    (
        const fvMesh&,
        const dictionary&
    );

    //- boundaries for named patches
    void namedBoundaries
    (
        const fvMesh&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&,
        bool strictPatchNameChecking
    );

    //- boundary dictionary entries for named patches
    void namedBoundaries
    (
        const fvMesh&,
        const dictionary&,
        bool strictPatchNameChecking
    );


public:

    // Constructors

        //- Construct dictionary
        boundarySetup
        (
            const fvMesh&,
            const word&,
            const dictionary&,
            bool strictPatchNameChecking = true,
            bool printBoundaries = false
        );

        //- Construct from components.
        boundarySetup
        (
            const fvMesh&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            bool strictPatchNameChecking = true,
            bool printBoundaries = false
        );

    // Destructor

        ~boundarySetup();


    // Member functions

    // Return dictionary
    const dictionary& boundaryDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "cfdTools/general/boundarySetup/boundarySetup.C"
#endif

#endif

// ************************************************************************* //
