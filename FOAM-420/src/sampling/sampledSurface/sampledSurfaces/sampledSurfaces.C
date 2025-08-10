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
    (c) 2016 OpenCFD Ltd.
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSurface/sampledSurfaces/sampledSurfaces.H"
#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "sampledSurface/sampledTriSurfaceMesh/sampledTriSurfaceMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurfaces, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSurfaces,
        dictionary
    );
}

bool Foam::sampledSurfaces::verbose_ = false;
Foam::scalar Foam::sampledSurfaces::mergeTol_ = 1e-10;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSurfaces::writeGeometry() const
{
    // Write to time directory under outputPath_
    // Skip surface without faces (eg, a failed cut-plane)


    fileName outputDir = outputPath_/time_.timeName();

    if (!writeTimeName_)
    {
        outputDir = outputPath_;
    }

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        if (Pstream::parRun())
        {
            if (Pstream::master() && mergedList_[surfI].size())
            {
                formatter_->write
                (
                    outputDir,
                    s.name(),
                    mergedList_[surfI]
                );
            }
        }
        else if (s.faces().size())
        {
            formatter_->write(outputDir, s.name(), s);
        }
    }
}


void Foam::sampledSurfaces::writeOriginalIds()
{
    const word fieldName = "Ids";
    const fileName outputDir = outputPath_/time_.timeName();

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        if (isA<sampledTriSurfaceMesh>(s))
        {
            const sampledTriSurfaceMesh& surf =
                dynamicCast<const sampledTriSurfaceMesh&>(s);

            if (surf.keepIds())
            {
                const labelList& idLst = surf.originalIds();

                Field<scalar> ids(idLst.size());
                forAll(idLst, i)
                {
                    ids[i] = idLst[i];
                }

                writeSurface(ids, surfI, fieldName, outputDir);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<sampledSurface>(),
    mesh_(refCast<const fvMesh>(meshObr_)),
    loadFromFiles_(false),
    outputPath_(fileName::null),
    fieldSelection_(),
    interpolationScheme_(word::null),
    writeTimeName_(dict.lookupOrDefault("writeTimePath", true)),
    writeStats_(dict.lookupOrDefault("writeStatistics", false)),
    writeFields_(dict.lookupOrDefault("writeFields", true)),
    mergedList_(),
    formatter_(nullptr)
{
    if (Pstream::parRun())
    {
        if (writeTimeName_)
        {
            outputPath_ = mesh_.time().path()/"../postProcessing";
        }
        else
        {
            outputPath_ = "../postProcessing";
        }
    }
    else
    {
        if (writeTimeName_)
        {
            outputPath_ = mesh_.time().path()/"postProcessing";
        }
        else
        {
            outputPath_ = "postProcessing";
        }
    }

    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

    outputPath_ = outputPath_/name;

    // Remove ".."
    outputPath_.clean();

    read(dict);

    if (!mesh_.conformal())
    {
        FatalIOErrorInFunction(dict)
            << "The " << type() << " function object is not compatible with "
            << "non-conformal coupling." << exit(FatalIOError);
    }
}


Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::regionFunctionObject(name, obr, dict),
    PtrList<sampledSurface>(),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    fieldSelection_(),
    interpolationScheme_(word::null),
    writeTimeName_(dict.lookupOrDefault("writeTimePath", true)),
    writeStats_(dict.lookupOrDefault("writeStatistics", false)),
    writeFields_(dict.lookupOrDefault("writeFields", true)),
    mergedList_(),
    formatter_(nullptr)
{

    read(dict);

    fileName dirName(dict.lookupOrDefault("dirname",fileName("postProcessing")));
    if (Pstream::parRun())
    {
        if (writeTimeName_)
        {
            outputPath_ = time_.path()/".."/dirName;
        }
        else
        {
            outputPath_ = ".."/dirName;
        }
    }
    else
    {
        if (writeTimeName_)
        {
            outputPath_ = time_.path()/dirName;
        }
        else
        {
            outputPath_ = dirName;
        }
    }

    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

    outputPath_ = outputPath_/name;

    // Remove ".."
    outputPath_.clean();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::~sampledSurfaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSurfaces::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::sampledSurfaces::execute()
{
    // Prevent surface to be updated twice
    // if write time and allow to update
    // when executeControl is set to true
    if (!mesh_.time().writeTime())
    {
        update();
    }

    return true;
}


bool Foam::sampledSurfaces::write()
{
    if (size())
    {
        // Finalize surfaces, merge points etc.
        update();

        const label nFields = classifyFields();

        if (Pstream::master())
        {
            if (debug)
            {
                if (writeTimeName_)
                {
                    Pout<< "Creating directory "
                        << outputPath_/time_.timeName() << nl << endl;
                }
                else
                {
                    Pout<< "Creating directory "
                        << outputPath_ << nl << endl;
                }
            }
            if (writeTimeName_)
            {
                if (!formatter_->sepFile())
                mkDir(outputPath_/time_.timeName());
            }
            else
            {
                if (!formatter_->sepFile())
                mkDir(outputPath_);
            }
        }
        // Write geometry first if required,
        // or when no fields would otherwise be written
        if (nFields == 0 || formatter_->separateGeometry())
        {
            writeGeometry();
        }

        const IOobjectList objects(obr(), obr().time().timeName());

        sampleAndWrite<volScalarField>(objects);
        sampleAndWrite<volVectorField>(objects);
        sampleAndWrite<volSphericalTensorField>(objects);
        sampleAndWrite<volSymmTensorField>(objects);
        sampleAndWrite<volTensorField>(objects);

        sampleAndWrite<surfaceScalarField>(objects);
        sampleAndWrite<surfaceVectorField>(objects);
        sampleAndWrite<surfaceSphericalTensorField>(objects);
        sampleAndWrite<surfaceSymmTensorField>(objects);
        sampleAndWrite<surfaceTensorField>(objects);
    }

    return true;
}


bool Foam::sampledSurfaces::read(const dictionary& dict)
{
    bool surfacesFound = dict.found("surfaces");

    if (surfacesFound)
    {
        dict.lookup("fields") >> fieldSelection_;

        dict.lookup("interpolationScheme") >> interpolationScheme_;
        const word writeType(dict.lookup("surfaceFormat"));

        // Define the surface formatter
        // Optionally defined extra controls for the output formats
        formatter_ = surfaceWriter::New
        (
            writeType,
            dict.subOrEmptyDict("formatOptions").subOrEmptyDict(writeType)
        );

        PtrList<sampledSurface> newList
        (
            dict.lookup("surfaces"),
            sampledSurface::iNew(mesh_)
        );
        transfer(newList);

        if (Pstream::parRun())
        {
            mergedList_.setSize(size());
        }

        // Ensure all surfaces and merge information are expired
        expire();

        if (this->size())
        {
            Info<< type() << " " << name() << ": "
                << "Reading surface description:" << nl;

            forAll(*this, surfI)
            {
                Info<< "    " << operator[](surfI).name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample surfaces:" << nl << "(" << nl;

        forAll(*this, surfI)
        {
            Pout<< "  " << operator[](surfI) << endl;
        }
        Pout<< ")" << endl;
    }

    return true;
}


void Foam::sampledSurfaces::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        expire();
    }

    // pointMesh and interpolation will have been reset in mesh.update
}


void Foam::sampledSurfaces::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        expire();
    }
}


void Foam::sampledSurfaces::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        expire();
    }
}


bool Foam::sampledSurfaces::needsUpdate() const
{
    forAll(*this, surfI)
    {
        if (operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::sampledSurfaces::expire()
{
    bool justExpired = false;

    forAll(*this, surfI)
    {
        if (operator[](surfI).expire())
        {
            justExpired = true;
        }

        // Clear merge information
        if (Pstream::parRun())
        {
            mergedList_[surfI].clear();
        }
    }

    // true if any surfaces just expired
    return justExpired;
}


bool Foam::sampledSurfaces::update()
{
    bool updated = false;

    if (!needsUpdate())
    {
        return updated;
    }

    // Serial: quick and easy, no merging required
    if (!Pstream::parRun())
    {
        forAll(*this, surfI)
        {
            if (operator[](surfI).update())
            {
                updated = true;
            }
        }

        return updated;
    }

    // Dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_*mesh_.bounds().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
            << mergeDim << " metre" << endl;
    }

    forAll(*this, surfI)
    {
        sampledSurface& s = operator[](surfI);

        if (s.update())
        {
            updated = true;
            mergedList_[surfI].merge(s, mergeDim);
        }
    }

    return updated;
}


Foam::scalar Foam::sampledSurfaces::mergeTol()
{
    return mergeTol_;
}


Foam::scalar Foam::sampledSurfaces::mergeTol(const scalar tol)
{
    scalar oldTol = mergeTol_;
    mergeTol_ = tol;
    return oldTol;
}


// ************************************************************************* //
