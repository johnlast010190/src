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
    (c) 2016-2020 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "ensight/mesh/ensightMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::options::options(IOstream::streamFormat format)
:
    format_(format),
    lazy_(false),
    internal_(true),
    boundary_(true),
    deprecatedOrder_(false),
    patchPatterns_(),
    faceZonePatterns_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOstream::streamFormat Foam::ensightMesh::options::format() const
{
    return format_;
}


bool Foam::ensightMesh::options::lazy() const
{
    return lazy_;
}


bool Foam::ensightMesh::options::useInternalMesh() const
{
    return internal_;
}


bool Foam::ensightMesh::options::useBoundaryMesh() const
{
    return boundary_;
}


bool Foam::ensightMesh::options::deprecatedOrder() const
{
    return deprecatedOrder_;
}


bool Foam::ensightMesh::options::useFaceZones() const
{
    return faceZonePatterns_.valid();
}


bool Foam::ensightMesh::options::usePatchSelection() const
{
    return boundary_ ? patchPatterns_.valid() : false;
}


void Foam::ensightMesh::options::reset()
{
    internal_ = true;
    boundary_ = true;
    patchPatterns_.clear();
    faceZonePatterns_.clear();
}


void Foam::ensightMesh::options::lazy(const bool b)
{
    lazy_ = b;
}

void Foam::ensightMesh::options::useInternalMesh(bool on)
{
    internal_ = on;
}

void Foam::ensightMesh::options::useBoundaryMesh(bool on)
{
    boundary_ = on;

    if (!boundary_)
    {
        if (patchPatterns_.valid())
        {
            patchPatterns_.clear();

            WarningInFunction
                << "Deactivating boundary, removed old patch selection"
                << endl;
        }
    }
}

void Foam::ensightMesh::options::deprecatedOrder(const bool b)
{
    deprecatedOrder_ = b;
}


void Foam::ensightMesh::options::patchSelection
(
    const wordReList& patterns
)
{
    if (!boundary_)
    {
        WarningInFunction
            << "Ignoring patch selection, boundary is disabled"
            << endl;
    }
    else
    {
        patchPatterns_.reset(new wordReList(patterns));
    }
}


void Foam::ensightMesh::options::faceZoneSelection
(
    const wordReList& patterns
)
{
    faceZonePatterns_.reset(new wordReList(patterns));
}


const Foam::wordReList& Foam::ensightMesh::options::patchSelection() const
{
    if (usePatchSelection())
    {
        return patchPatterns_();
    }
    else
    {
        return wordReList::null();
    }
}


const Foam::wordReList& Foam::ensightMesh::options::faceZoneSelection() const
{
    if (faceZonePatterns_.valid())
    {
        return faceZonePatterns_();
    }
    else
    {
        return wordReList::null();
    }
}


// ************************************************************************* //
