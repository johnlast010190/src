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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"

#include "db/Time/Time.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "surfMesh/surfMesh.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/meshShapes/traits/faceTraits.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::MeshedSurfaceProxy<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


template<class Face>
bool Foam::MeshedSurfaceProxy<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes(), ext, verbose, "writing"
    );
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const fileName& name,
    const MeshedSurfaceProxy& surf
)
{
    write(name, name.ext(), surf);
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const fileName& name,
    const word& ext,
    const MeshedSurfaceProxy& surf
)
{
    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }

    auto mfIter = writefileExtensionMemberFunctionTablePtr_->find(ext);

    if (!mfIter.found())
    {
        FatalErrorInFunction
            << "Unknown file extension " << ext << nl << nl
            << "Valid types:" << nl
            << flatOutput(writeTypes().sortedToc()) << nl
            << exit(FatalError);
    }

    mfIter()(name, surf);
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    // the surface name to be used
    const word name(surfName.size() ? surfName : surfaceRegistry::defaultName);

    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }


    // The local location
    const fileName objectDir
    (
        t.timePath()/surfaceRegistry::prefix/name/surfMesh::meshSubDir
    );

    if (!isDir(objectDir))
    {
        mkDir(objectDir);
    }


    // write surfMesh/points
    {
        pointIOField io
        (
            IOobject
            (
                "points",
                t.timeName(),
                surfMesh::meshSubDir,
                t,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        OFstream os
        (
            objectDir/io.name(),
            t.writeFormat(),
            IOstream::currentVersion,
            t.writeCompression()
        );

        io.writeHeader(os);

        os  << this->points();

        io.writeEndDivider(os);
    }


    // write surfMesh/faces
    {
        faceCompactIOList io
        (
            IOobject
            (
                "faces",
                t.timeName(),
                surfMesh::meshSubDir,
                t,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        OFstream os
        (
            objectDir/io.name(),
            t.writeFormat(),
            IOstream::currentVersion,
            t.writeCompression()
        );
        io.writeHeader(os);

        if (this->useFaceMap())
        {
            // this is really a bit annoying (and wasteful) but no other way
            os  << reorder(this->faceMap(), this->surfFaces());
        }
        else
        {
            os  << this->surfFaces();
        }

        io.writeEndDivider(os);
    }


    // write surfMesh/surfZones
    {
        surfZoneIOList io
        (
            IOobject
            (
                "surfZones",
                t.timeName(),
                surfMesh::meshSubDir,
                t,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // write as ascii
        OFstream os(objectDir/io.name());
        io.writeHeader(os);

        os  << this->surfZones();

        io.writeEndDivider(os);
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::MeshedSurfaceProxy
(
    const pointField& pointLst,
    const List<Face>& faceLst,
    const List<surfZone>& zoneLst,
    const List<label>& faceMap
)
:
    points_(pointLst),
    faces_(faceLst),
    zones_(zoneLst),
    faceMap_(faceMap)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::~MeshedSurfaceProxy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
inline Foam::label Foam::MeshedSurfaceProxy<Face>::nTriangles() const
{
    if (faceTraits<Face>::isTri())
    {
        return this->size();
    }
    else
    {
        label nTri = 0;
        const List<Face>& faceLst = this->surfFaces();
        forAll(faceLst, facei)
        {
            nTri += faceLst[facei].nTriangles();
        }

        return nTri;
    }
}


// ************************************************************************* //
