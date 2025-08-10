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

InNamespace
    Foam

Description
    Bound the given scalar field if it has gone unbounded.

    Used extensively in RAS and LES turbulence models, but also of use
    within solvers.

SourceFiles
    bound.C

\*---------------------------------------------------------------------------*/

#ifndef bound_H
#define bound_H

#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Bound the given scalar field if it has gone unbounded.
//  Return the bounded field.
//  Used extensively in RAS and LES turbulence models.
volScalarField& bound
(
    volScalarField&,
    const dimensionedScalar& lowerBound,
    bool warn = true
);

//- Bound the given scalar field if it has gone unbounded.
//  Return the bounded field.
//  Used extensively in RAS and LES turbulence models.
volScalarField& bound
(
    volScalarField&,
    const dimensionedScalar& lowerBound,
    const dimensionedScalar& upperBound,
    bool warn = true
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
