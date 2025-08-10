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

Namespace
    Foam::blockMeshTools

Description
    Tools for parsing label(List) with dictionary lookup.

SourceFiles
    blockMeshTools.C
    blockMeshToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef blockMeshTools_H
#define blockMeshTools_H

#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockMeshTools
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- In-place read with dictionary lookup
    void read(Istream&, label&, const dictionary&);

    //- In-place read with dictionary lookup
    template<class T>
    void read(Istream&, List<T>&, const dictionary&);

    //- Return-read with dictionary lookup
    label read(Istream&, const dictionary&);

    //- Return-read with dictionary lookup
    template<class T>
    List<T> read(Istream& is, const dictionary&);

    //- Write with dictionary lookup
    void write(Ostream&, const label, const dictionary&);

    //- Linear search for label entry
    const keyType& findEntry(const dictionary&, const label);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockMeshTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blockMeshTools/blockMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
