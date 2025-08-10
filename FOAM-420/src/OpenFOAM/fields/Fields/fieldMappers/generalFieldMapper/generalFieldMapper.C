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
    (c) 2019-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/fieldMappers/generalFieldMapper/generalFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::generalFieldMapper::directAddressing() const
{
    FatalErrorInFunction
        << "attempt to access null direct addressing"
        << abort(FatalError);

    return labelUList::null();
}


const Foam::labelUList& Foam::generalFieldMapper::indirectAddressing() const
{
    FatalErrorInFunction
        << "attempt to access null indirect addressing"
        << abort(FatalError);

    return labelUList::null();
}


const Foam::labelListList& Foam::generalFieldMapper::addressing() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation addressing"
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::generalFieldMapper::weights() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation weights"
        << abort(FatalError);

    return scalarListList::null();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, , generalFieldMapper)


DEFINE_FIELD_MAPPER_OPERATOR(label, , generalFieldMapper)


forAllVectorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, generalFieldMapper)

forAllTensorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, generalFieldMapper)

forAllDiagTensorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, generalFieldMapper)

forAllSphericalTensorNTypes(DEFINE_FIELD_MAPPER_OPERATOR, generalFieldMapper)


// ************************************************************************* //
