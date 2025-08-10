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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "volumeHeight.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(volumeHeight, 0);
    addToRunTimeSelectionTable(fieldInit, volumeHeight, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

Foam::scalar Foam::fieldInitializations::volumeHeight::compVolume
(
    const labelList& cellIds,
    const List<scalar>& heights,
    const scalar maxHeight
) const
{
    scalar vloc(0.0);

    forAll(cellIds, ci)
    {
        if (heights[ci] <= maxHeight)
        {
            vloc += mesh().V()[cellIds[ci]];
        }
    }

    return returnReduce(vloc, sumOp<scalar>());
}


void Foam::fieldInitializations::volumeHeight::setFieldValue
(
    const labelList& cellIds,
    const List<scalar>& heights,
    const scalar minHeight,
    const scalar maxHeight,
    const scalar setValue
) const
{
    volScalarField& f
        = const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()));

    forAll(cellIds, ci)
    {
        if
        (
            heights[ci] <= maxHeight
         && heights[ci] > minHeight
        )
        {
            f[cellIds[ci]] = setValue;
        }
    }
}


bool Foam::fieldInitializations::volumeHeight::setCellValue
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

        uniformDimensionedVectorField g
        (
            IOobject
            (
                "g",
                mesh().time().constant(),
                localDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

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

            // get size of list over all processors
            const label nCells = returnReduce(ccells.toc().size(), sumOp<label>());

            // read target volume/height
            scalar vtarget(0.0);
            scalar htarget(-GREAT);
            if (setDict.found("fillVolume"))
            {
                vtarget = readScalar(setDict.lookup("fillVolume"));
                if (vtarget < VSMALL)
                {
                    FatalErrorInFunction
                        << "Specified fillVolume must be greater than zero. "
                        << exit(FatalError);
                }
            }
            else if (setDict.found("fillHeight"))
            {
                htarget = readScalar(setDict.lookup("fillHeight"));
            }
            else
            {
                WarningInFunction
                    << "No fillVolume or fillHeight specified, doing nothing"
                    << nl << endl;

                continue;
            }

            // create lists for cellid and height
            labelList cells = ccells.toc();
            List<scalar> heights(cells.size());

            // make sure cell centres are initialised
            const volVectorField& C = mesh().C();

            // compute heights
            forAll(cells, ci)
            {
                label celli = cells[ci];
                heights[ci] = -C[celli] & g[celli]/mag(g[celli]);
            }

            // get global min/max height
            scalar minHeight(min(heights));
            scalar maxHeight(max(heights));
            reduce(
                std::tie(minHeight, maxHeight),
                ParallelOp<minOp<scalar>, maxOp<scalar>>{}
            );

            if (htarget != -GREAT)
            {
                htarget += minHeight;
                setFieldValue(cells,heights,-GREAT,htarget);
                continue;
            }

            // get total volume
            scalar vtot = compVolume(cells,heights,GREAT);

            // if target volume > cellSet volume
            // fill everything and throw warning
            if (vtarget >= vtot)
            {
                setFieldValue(cells,heights);

                WarningInFunction
                    << "Specified target fill volume " << vtarget
                    << " is larger than the total cellSet volume "
                    << vtot << nl << endl;

                continue;
            }
            // get initial height estimate
            scalar h(minHeight + (maxHeight - minHeight)/(vtot - vtarget));

            // search algorithm to find height
            bool found(false);
            scalar hold(minHeight);
            scalar holdold(minHeight);
            scalar vcur = compVolume(cells,heights,h);
            scalar vold(vcur);

            // search tolerance
            scalar mindiffdef = cbrt(vtot/nCells)/5.0;
            scalar mindiff = setDict.lookupOrDefault<scalar>("minDiff", mindiffdef);

            while (!found)
            {
                scalar hdiffmin = mindiff;

                if (mag(h-hold) > hdiffmin)
                {
                    // do bisection search
                    if (vcur < vtarget)
                    {
                        minHeight = h;
                    }
                    else
                    {
                        maxHeight = h;
                    }

                    holdold = hold;
                    hold = h;
                    h = (minHeight+maxHeight)/2.0;
                }
                else
                {
                    // do linear search
                    if (vcur < vtarget)
                    {
                        hdiffmin = mindiff;
                    }
                    else
                    {
                        hdiffmin = -mindiff;
                    }

                    holdold = hold;
                    hold = h;
                    h += hdiffmin;
                }

                vold = vcur;
                vcur = compVolume(cells,heights,h);

                // search end condition
                if
                (
                    h == holdold
                 && min(vcur,vold) <= vtarget
                 && max(vcur,vold) >= vtarget
                )
                {
                    found = true;

                    // compute volume of completely filled cells and set value
                    scalar vfill = compVolume(cells,heights,min(h, hold));
                    setFieldValue(cells,heights,-GREAT,min(h, hold));

                    // compute partial fill value from missing volume and set value
                    scalar alphadiff = (vtarget - vfill)/max(max(vcur,vold) - vfill, SMALL);
                    setFieldValue(cells,heights,min(h, hold),max(h, hold),alphadiff);

                    if (debug)
                    {
                        Info<< "height : " << h
                             << nl << "vcur : " << vcur
                             << nl << "vfill : " << vfill
                             << nl << "mindiff : " << mindiff
                             << nl << "alphadiff : " << alphadiff << endl;
                    }
                }
            }

            // check and report total phase volume
            scalar checkVolume(0.0);
            forAll(cells, ci)
            {
                label celli = cells[ci];
                checkVolume += f[celli]*mesh().V()[celli];
            }
            reduce(checkVolume, sumOp<scalar>());
            Info<< "total filled volume " << checkVolume << endl;
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

Foam::fieldInitializations::volumeHeight::volumeHeight
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

void Foam::fieldInitializations::volumeHeight::correct()
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

    setCellValue(setDicts, cellSets);

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
