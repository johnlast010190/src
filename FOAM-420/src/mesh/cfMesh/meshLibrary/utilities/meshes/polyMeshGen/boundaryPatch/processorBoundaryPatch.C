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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGen/boundaryPatch/processorBoundaryPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/IOstreams/Istream.H"
#include "db/IOstreams/token/token.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorBoundaryPatch, 0);
addToRunTimeSelectionTable(boundaryPatchBase, processorBoundaryPatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

processorBoundaryPatch::processorBoundaryPatch
(
    const word& name,
    const word& type,
    const label nFaces,
    const label startFace,
    const label myProcNo,
    const label neighbProcNo
)
:
    boundaryPatchBase(name, type, nFaces, startFace),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo)
{}

processorBoundaryPatch::processorBoundaryPatch
(
    const word& name,
    const dictionary& dict
)
:
    boundaryPatchBase(name, dict),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo")))
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary processorBoundaryPatch::dict() const
{
    dictionary dict;

    dict.add("type", type_);

    dict.add("nFaces", nFaces_);
    dict.add("startFace", startFace_);
    dict.add("myProcNo", myProcNo_);
    dict.add("neighbProcNo", neighbProcNo_);

    return dict;
}

void processorBoundaryPatch::write(Ostream& os) const
{
    this->operator<<(os);
}

void processorBoundaryPatch::writeDict(Ostream& /*os*/) const
{

}

Ostream& processorBoundaryPatch::operator<<(Ostream& os) const
{
    os  << patchName() << nl << token::BEGIN_BLOCK << nl
        << "    type         " << patchType() << token::END_STATEMENT << nl
        << "    nFaces       " << patchSize() << token::END_STATEMENT << nl
        << "    startFace    " << patchStart() << token::END_STATEMENT << nl
        << "    myProcNo     " << myProcNo_ << token::END_STATEMENT << nl
        << "    neighbProcNo " << neighbProcNo_
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;

    return os;
}

Istream& processorBoundaryPatch::operator>>(Istream& is)
{
    token t;
    is >> name_ >> t;
    is >> t >> type_ >> t;
    is >> t >> nFaces_ >> t;
    is >> t >> startFace_ >> t;
    is >> t >> myProcNo_ >> t;
    is >> t >> neighbProcNo_ >> t;
    is >> t;

    return is;
}

bool processorBoundaryPatch::operator!=(const processorBoundaryPatch& wp) const
{
    if (name_ != wp.name_)
    {
        return true;
    }
    else if (type_ != wp.name_)
    {
        return true;
    }
    else if ((nFaces_ != wp.nFaces_) || (startFace_ != wp.startFace_))
    {
        return true;
    }
    else if (
        (myProcNo_ != wp.myProcNo_) || (neighbProcNo_ != wp.neighbProcNo_)
    )
    {
        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
