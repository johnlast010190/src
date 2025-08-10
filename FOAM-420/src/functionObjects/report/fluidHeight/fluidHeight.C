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
    (c) 2011 OpenFOAM Foundation
    (c) 2017-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fluidHeight/fluidHeight.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fluidHeight, 0);
    addToRunTimeSelectionTable(functionObject, fluidHeight, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fluidHeight::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Region fluid heights");
    writeCommented(os, "Time");
    forAll(reportFaces_, regioni)
    {
        writeDelimited(os, regionNames_[regioni]);
    }
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fluidHeight::fluidHeight
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    vofName_(dict.lookup("fieldName")),
    tValue_(dict.lookupOrDefault<scalar>("threshold", 0.5)),
    reportCells_((refCast<const polyMesh>(obr_)).nCells(), 0),
    reportFaces_(),
    regionNames_()
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fluidHeight::~fluidHeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fluidHeight::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    vofName_ = word(dict.lookup("fieldName"));
    tValue_ = dict.lookupOrDefault<scalar>("threshold", 0.5);

    //select cells
    //if sets entry is not present, select all cells
    if (dict.found("sets"))
    {
        PtrList<entry> regions(dict.lookup("sets"));

        // number of regions
        label nRegions = regions.size();
        reportFaces_.resize(nRegions);
        regionNames_.resize(nRegions);

        forAll(regions, regioni)
        {
            const entry& region = regions[regioni];

            autoPtr<topoSetSource> cellSelector =
                topoSetSource::New(region.keyword(), mesh_, region.dict());

            cellSet selectedCellSet
            (
                mesh_,
                "cellSet",
                mesh_.nCells()/10+1  // Reasonable size estimate.
            );

            cellSelector->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            regionNames_[regioni] = word(region.dict().lookup("name"));

            // add faces to list
            forAllConstIter(labelHashSet, selectedCellSet, iter)
            {
                reportCells_.set(iter.key(), 1);
                const cell& cFaces = mesh_.cells()[iter.key()];
                forAll(cFaces, cFaceI)
                {
                    reportFaces_[regioni].set(cFaces[cFaceI], 1);
                }
            }
        }
    }
    else // add all cells/faces to list
    {
        forAll(mesh_.C(), cI)
        {
            reportCells_.set(cI, 1);
        }

        reportFaces_.resize(1);
        regionNames_.resize(1);
        regionNames_[0] = word("region0");
        forAll(mesh_.C(), cellI)
        {
            const cell& cFaces = mesh_.cells()[cellI];
            forAll(cFaces, cFaceI)
            {
                reportFaces_[0].set(cFaces[cFaceI], 1);
            }
        }
    }
    return true;
}


bool Foam::functionObjects::fluidHeight::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;
    writeTime(file());

    // compute gravitation direction vector
    const uniformDimensionedVectorField& g =
        obr_.lookupObject<uniformDimensionedVectorField>("g");
    vector upDir = -g.value() / mag(g.value());

    // get reference to two-phase volume fraction field
    const volScalarField& alpha = obr_.lookupObject<volScalarField>(vofName_);

    const volVectorField& CC(mesh_.C());

    // store cell center values on face and sort for min/max
    string localMin("localMin");
    surfaceScalarField faceMin( fvc::interpolate(alpha, IStringStream(localMin)()) );
    string localMax("localMax");
    surfaceScalarField faceMax( fvc::interpolate(alpha, IStringStream(localMax)()) );

    // get boundary field data
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    tmp<surfaceVectorField> tdelta(mesh_.delta());
    const surfaceVectorField& delta = tdelta();

    //loop through all marked cells
    forAll(reportFaces_, regioni)
    {
        // init fluid height
        scalar maxHeight = -GREAT;

        forAll(reportFaces_[regioni], fI)
        {
            if (reportFaces_[regioni].get(fI) == 1)
            {
                // check owner and neighbour for boundary
                if (fI < mesh_.nInternalFaces())
                {
                    if (faceMax[fI] >= tValue_ && faceMin[fI] < tValue_)
                    {
                        scalar iFrac
                            = (faceMax[fI] - tValue_)
                            /(faceMax[fI] - faceMin[fI]);

                        label maxCellI = -1;
                        label minCellI = -1;

                        if (alpha[mesh_.owner()[fI]] == faceMax[fI])
                        {
                            maxCellI = mesh_.owner()[fI];
                            minCellI = mesh_.neighbour()[fI];
                        }
                        else
                        {
                            minCellI = mesh_.owner()[fI];
                            maxCellI = mesh_.neighbour()[fI];
                        }
                        vector maxC = CC[maxCellI];
                        vector max2minC = CC[minCellI] - maxC;
                        vector interfacePos = maxC + iFrac*max2minC;
                        maxHeight = max(maxHeight, upDir & interfacePos);
                    }
                }
                else // processor patch?
                {
                    label patchI = patches.whichPatch(fI);
                    if
                    (
                        isA<directPolyPatch>(patches[patchI])
                      && alpha.boundaryField()[patchI].size()
                    )
                    {
                        label pfI = patches[patchI].whichFace(fI);
                        if
                        (
                            !patches[patchI].coupled()
                        )
                        {
                            label cellID = mesh_.faceOwner()[fI];
                            if (alpha[cellID] >= tValue_)
                            {
                                maxHeight = max
                                (
                                    maxHeight,
                                    (
                                        upDir
                                      & mesh_.Cf().boundaryField()[patchI][pfI]
                                    )
                                );
                            }
                        }
                        else
                        {
                            scalar fMax = faceMax.boundaryField()[patchI][pfI];
                            scalar fMin = faceMin.boundaryField()[patchI][pfI];
                            if
                            (
                                fMax >= tValue_
                                && fMin < tValue_
                            )
                            {
                                label cI =  mesh_.faceOwner()[fI];
                                scalar iFrac = (fMax - tValue_)/(fMax - fMin);

                                if (alpha[cI] == fMax)
                                {
                                    vector interfacePos =
                                        CC[cI] + iFrac*delta[fI];

                                    maxHeight =
                                        max(maxHeight,upDir & interfacePos);
                                }
                            }
                        }
                    }
                }
            }
        }

        reduce(maxHeight, maxOp<scalar>());

        Log << tab << "Maximum fluid height in region "
             << regionNames_[regioni] << ": " << maxHeight << " [m]." << endl;
        file() << token::TAB << maxHeight;
    }
    file()<< endl;


    return true;
}


bool Foam::functionObjects::fluidHeight::write()
{
    return true;
}


// ************************************************************************* //
