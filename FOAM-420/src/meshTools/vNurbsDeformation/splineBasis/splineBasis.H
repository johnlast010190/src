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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    splineBasis

Description
    B-spline basis function

\*---------------------------------------------------------------------------*/

#ifndef splineBasis_H
#define splineBasis_H

#include "db/dictionary/dictionary.H"
#include "fields/Fields/scalarField/scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class splineBasis Declaration
\*---------------------------------------------------------------------------*/

class splineBasis
{
    // private data

    // interpolation points, degree, knots
    label nCPs, nDegree, nKnots;

    // knot array
    scalarField knots;

    // private member functions
    void createKnots();

public:

    // Constructor

    splineBasis
    (
        const label nPoints,
        const label deg
    );
    // copy constructor
    splineBasis
    (
        const splineBasis& basis
    );

    ~splineBasis()
    {}

    // public member functions

    // basis value
    scalar value
    (
        const label& iCP,
        const label& deg,
        const scalar& u
    ) const;

    // basis der
    scalar derivative
    (
        const label& iCP,
        const label& deg,
        const scalar& u
    ) const;

    // find knot span
    label span(const scalar& u) const;

    // get knots
    const scalarField& getKnots() const;

    void changeKnots
    (
        const scalarField& k
    );

    // get private data
    const label& nbPoints() const;
    const label& nbKnots() const;
    const label& degree() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
