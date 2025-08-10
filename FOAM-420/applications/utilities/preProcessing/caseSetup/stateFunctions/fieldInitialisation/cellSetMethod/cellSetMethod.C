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
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(cellSetMethod, 0);
    addToRunTimeSelectionTable(fieldInit, cellSetMethod, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

template<>
bool Foam::fieldInitializations::cellSetMethod::setCellValue<Foam::scalar>
(
    const PtrList<entry>& setDicts,
    const PtrList<cellSet>& cellSets
) const
{

    if (localDb().foundObject<volScalarField>(name()))
    {
        volScalarField& f
            = const_cast<volScalarField&>
            (localDb().lookupObject<volScalarField>(name()));

        if (initDict().found("defaultValue"))
        {
            f.primitiveFieldRef() = Field<scalar>
            (
                word("defaultValue"), initDict(), f.size()
            );
        }

        forAll(setDicts, seti)
        {
            const dictionary& setDict = setDicts[seti].dict();
            const cellSet& ccells = cellSets[seti];

            scalar v = readScalar(setDict.lookup("value"));

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::cellSetMethod::cellSetMethod
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::cellSetMethod::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    PtrList<entry> setDicts(initDict().lookup("setSources"));

    PtrList<cellSet> cellSets(setDicts.size());

    forAll(setDicts, seti)
    {
        const entry& setDict = setDicts[seti];

        autoPtr<topoSetSource> cellSelector =
            topoSetSource::New(setDict.keyword(), mesh(), setDict.dict());

        cellSets.set
        (
            seti,
            new cellSet
            (
                mesh(),
                word("cellSet") + Foam::name(seti),
                mesh().nCells()
            )
        );

        cellSelector->applyToSet
        (
            topoSetSource::ADD,
            cellSets[seti]
        );
    }


    if (setCellValue<scalar>(setDicts, cellSets));
    else if (setCellValue<vector>(setDicts, cellSets));
    else if (setCellValue<tensor>(setDicts, cellSets));
    else if (setCellValue<symmTensor>(setDicts, cellSets));
    else if (setCellValue<sphericalTensor>(setDicts, cellSets));
    else
    {
        WarningInFunction
            << "Could not find field " << name()
            << " in database. Initialisation cancelled."
            << endl;
    }

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
