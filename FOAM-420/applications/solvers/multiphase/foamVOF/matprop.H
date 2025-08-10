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
    (c) 2019 Esi Ltd.

Application
    foamVOF

Group
    grpMultiphaseSolvers

Description
    Solver for multi-phase flows with two or more fluids. This file contains
    implementation details of the member functions of the MatProp class

\*---------------------------------------------------------------------------*/

#include "matprop.h"
using namespace Foam;

MatProp::MatProp(const Time& runTime, const dynamicFvMesh& mesh)
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
   // MatProp prop(transportProperties);
    cp1_=scalar(transportProperties.lookupOrDefault<scalar>("cp1", 6200));
    cp2_=scalar(transportProperties.lookupOrDefault<scalar>("cp2", 1007));
    k1_=scalar(transportProperties.lookupOrDefault<scalar>("k1", 0.6));
    k2_=scalar(transportProperties.lookupOrDefault<scalar>("k2", 0.025));
    rho1_=scalar(transportProperties.lookupOrDefault<scalar>("rho1", 1000));
    rho2_=scalar(transportProperties.lookupOrDefault<scalar>("rho2", 1.19));
    Sch_=scalar(transportProperties.lookupOrDefault<scalar>("Schmidt", 0.7));
    Pr_=scalar(transportProperties.lookupOrDefault<scalar>("Prandtl", 10));
    psat_=scalar(transportProperties.lookupOrDefault<scalar>("psat", 2300));

}

MatProp::MatProp(dictionary &transportProperties)
:
cp1_(transportProperties.lookupOrDefault<scalar>("cp1", 6200)),
cp2_(transportProperties.lookupOrDefault<scalar>("cp2", 1007)),
k1_(transportProperties.lookupOrDefault<scalar>("k1", 0.6)),
k2_(transportProperties.lookupOrDefault<scalar>("k2", 0.025))
{}



