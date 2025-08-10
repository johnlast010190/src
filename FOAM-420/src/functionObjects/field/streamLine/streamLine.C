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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "streamLine/streamLine.H"
#include "streamLine/streamLineParticleCloud.H"
#include "sampledSet/sampledSet/sampledSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(streamLine, 0);
    addToRunTimeSelectionTable(functionObject, streamLine, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::streamLine::track()
{
    IDLList<streamLineParticle> initialParticles;
    streamLineParticleCloud particles
    (
        mesh_,
        cloudName_,
        initialParticles
    );

    const sampledSet& seedPoints = sampledSetPtr_();

    forAll(seedPoints, i)
    {
        particles.addParticle
        (
            new streamLineParticle
            (
                mesh_,
                seedPoints[i],
                seedPoints.cells()[i],
                lifeTime_
            )
        );
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    Log << "    seeded " << nSeeds << " particles" << endl;

    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolation<scalar>> vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolation<vector>> vvInterp;

    label UIndex = -1;

    initInterpolations
    (
        nSeeds,
        UIndex,
        vsFlds,
        vsInterp,
        vvFlds,
        vvInterp
    );

    // Additional particle info
    streamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction?
        nSubCycle_,     // automatic track control:step through cells in steps?
        trackLength_,   // fixed track length

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::streamLine::streamLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    streamLineBase(name, runTime, dict)
{
    read(dict_);

    if (!mesh_.conformal())
    {
        FatalIOErrorInFunction(dict)
            << "The " << type() << " function object is not compatible with "
            << "non-conformal coupling." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::streamLine::~streamLine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::streamLine::read(const dictionary& dict)
{
    if (streamLineBase::read(dict))
    {
        bool subCycling = dict.found("nSubCycle");
        bool fixedLength = dict.found("trackLength");

        if (subCycling && fixedLength)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot both specify automatic time stepping (through '"
                << "nSubCycle' specification) and fixed track length (through '"
                << "trackLength')"
                << exit(FatalIOError);
        }

        nSubCycle_ = 1;
        if (dict.readIfPresent("nSubCycle", nSubCycle_))
        {
            trackLength_ = VGREAT;
            nSubCycle_ = max(nSubCycle_, 1);


            Info<< "    automatic track length specified through"
                << " number of sub cycles : " << nSubCycle_ << nl
                << endl;
        }
    }
    return true;
}


// ************************************************************************* //
