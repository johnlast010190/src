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
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cellSetMethod.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldInitializations::cellSetMethod::setCellValue
(
    const PtrList<entry>& setDicts,
    const PtrList<cellSet>& cellSets
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    if (localDb().foundObject<GeoField>(name()))
    {
        GeoField& f
            = const_cast<GeoField&>
            (localDb().lookupObject<GeoField>(name()));

        if (initDict().found("defaultValue"))
        {
            f.primitiveFieldRef() = Field<Type>
            (
                word("defaultValue"), initDict(), f.size()
            );
        }

        forAll(setDicts, seti)
        {
            const dictionary& setDict = setDicts[seti].dict();
            const cellSet& ccells = cellSets[seti];

            Type v = setDict.lookup("value");

            labelList cells = ccells.toc();
            forAll(cells, ci)
            {
                label celli = cells[ci];
                f[celli] = v;
            }
        }

        // set patches to adjacent field values
        forAll(f.boundaryField(), patchi)
        {
            f.boundaryFieldRef()[patchi] =
                f.boundaryField()[patchi].patchInternalField();
        }

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
