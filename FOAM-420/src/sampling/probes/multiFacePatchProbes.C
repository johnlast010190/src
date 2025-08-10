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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "probes/multiFacePatchProbes.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiFacePatchProbes, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiFacePatchProbes::findElements(const fvMesh& mesh)
{
    //zero probe based patch faces and areas
    probeFaces_.setSize(this->size());
    forAll(probeFaces_, probeI)
    {
        probeFaces_[probeI].setSize(0);
    }

    probeArea_.setSize(this->size());
    probeArea_ = 0;

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    label patchI = bm.findPatchID(patchName_);

    if (patchI == -1)
    {
        FatalErrorInFunction
            << " Unknown patch name "
            << patchName_ << endl
            << exit(FatalError);
    }

    const polyPatch& pp = bm[patchI];

    if (pp.size() > 0)
    {
        forAll(probeLocations(), probeI)
        {
            const point sample = probeLocations()[probeI];

            scalarField distSqr(magSqr(pp.faceCentres() - sample));

            DynamicList<label> cProbeFaces(10);

            forAll(distSqr, pi)
            {
                if (distSqr[pi] < radius2_)
                {
                    cProbeFaces.append(pi);
                    probeArea_[probeI] += pp.magFaceAreas()[pi];
                }
            }

            cProbeFaces.shrink();
            probeFaces_[probeI] = cProbeFaces;
        }
    }

    Pstream::listCombineReduce(probeArea_, plusOp<scalar>());

    forAll(probeArea_, paI)
    {
        if (probeArea_[paI] == 0)
        {
            FatalErrorInFunction
                << this->type() << " functionObject:" << name() << nl
                << "Zero face area detected for probe "
                << paI << " on patch " << patchName_ << " (" << patchI << ")"
                << " at location " << this->operator[](paI)
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiFacePatchProbes::multiFacePatchProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    probes(name, runTime, dict, loadFromFiles),
    patchName_(),
    radius2_(0),
    probeScale_(1),
    probeFaces_(),
    probeArea_()
{
    // When constructing probes above it will have called the
    // probes::findElements (since the virtual mechanism not yet operating).
    // Not easy to workaround (apart from feeding through flag into constructor)
    // so clear out any cells found for now.
    elementList_.clear();
    faceList_.clear();
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiFacePatchProbes::~multiFacePatchProbes()
{}


bool Foam::multiFacePatchProbes::write()
{
    if (this->size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);

        if (debug)
        {
            forAll(probeArea_, probeI)
            {
                Info<< "xyz: " << this->operator[](probeI)
                     << ", area: " << probeArea_[probeI] << endl;
            }
        }
    }

    return true;
}

bool Foam::multiFacePatchProbes::read(const dictionary& dict)
{
    dict.lookup("patchName") >> patchName_;
    radius2_ = sqr(readScalar(dict.lookup("probeRadius")));
    probeScale_ = dict.lookupOrDefault<scalar>("probeScale", 1);

    dict.lookup("probeLocations") >> *this;
    this->operator*=(probeScale_);
    dict.lookup("fields") >> fieldSelection_;
    findElements(mesh_);
    probes::prepare();

    return true;
}


// ************************************************************************* //
