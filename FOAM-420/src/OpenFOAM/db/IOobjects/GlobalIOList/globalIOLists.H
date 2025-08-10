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
    (c) 2015-17 OpenFOAM Foundation

Typedef
    Foam::globalIOLists

Description
    Typedefs for globalIOLists of primitive types. These are fully compatible
    with 'normal' IOLists except have global filePath() scope.

\*---------------------------------------------------------------------------*/

#ifndef globalIOLists_H
#define globalIOLists_H

#include "fields/Fields/primitiveFields.H"
#include "db/IOobjects/GlobalIOList/GlobalIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef GlobalIOList<label> labelGlobalIOList;
    typedef GlobalIOList<scalar> scalarGlobalIOList;
    typedef GlobalIOList<vector> vectorGlobalIOList;
    typedef GlobalIOList<sphericalTensor> sphericalTensorGlobalIOList;
    typedef GlobalIOList<symmTensor> symmTensorGlobalIOList;
    typedef GlobalIOList<tensor> tensorGlobalIOList;

    //- Template function for obtaining global status
    template<>
    inline bool typeGlobal<labelGlobalIOList>()
    {
        return true;
    }

    template<>
    inline bool typeGlobal<scalarGlobalIOList>()
    {
        return true;
    }

    template<>
    inline bool typeGlobal<vectorGlobalIOList>()
    {
        return true;
    }

    template<>
    inline bool typeGlobal<sphericalTensorGlobalIOList>()
    {
        return true;
    }

    template<>
    inline bool typeGlobal<symmTensorGlobalIOList>()
    {
        return true;
    }

    template<>
    inline bool typeGlobal<tensorGlobalIOList>()
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
