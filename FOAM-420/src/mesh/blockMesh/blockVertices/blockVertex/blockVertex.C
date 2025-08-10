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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "blockVertices/blockVertex/blockVertex.H"
#include "blockVertices/pointVertex/pointVertex.H"
#include "blockMeshTools/blockMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockVertex, 0);
    defineRunTimeSelectionTable(blockVertex, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockVertex::blockVertex()
{}


Foam::autoPtr<Foam::blockVertex> Foam::blockVertex::clone() const
{
    NotImplemented;
}


Foam::autoPtr<Foam::blockVertex> Foam::blockVertex::New
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction << "Constructing blockVertex" << endl;
    }

    token firstToken(is);

    if (firstToken.isPunctuation() && firstToken.pToken() == token::BEGIN_LIST)
    {
        // Putback the opening bracket
        is.putBack(firstToken);

        return autoPtr<blockVertex>
        (
            new blockVertices::pointVertex(dict, index, geometry, is)
        );
    }
    else if (firstToken.isWord())
    {
        const word faceType(firstToken.wordToken());

        const auto ctor = ctorTableLookup("blockVertex type", IstreamConstructorTable_(), faceType);
        return autoPtr<blockVertex>(ctor(dict, index, geometry, is));
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <word> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);

        return autoPtr<blockVertex>(nullptr);
    }
}


Foam::label Foam::blockVertex::read(Istream& is, const dictionary& dict)
{
    const dictionary* varDictPtr = dict.subDictPtr("namedVertices");
    if (varDictPtr)
    {
        return blockMeshTools::read(is, *varDictPtr);
    }
    return readLabel(is);
}


void Foam::blockVertex::write
(
    Ostream& os,
    const label val,
    const dictionary& d
)
{
    const dictionary* varDictPtr = d.subDictPtr("namedVertices");
    if (varDictPtr)
    {
        blockMeshTools::write(os, val, *varDictPtr);
    }
    else
    {
        os << val;
    }
}


// ************************************************************************* //
