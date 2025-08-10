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
    (c) 2018-2021 Esi Ltd.
    (c) 2012-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "solidRotationSource.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvm/fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidRotationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        solidRotationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidRotationSource::solidRotationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    GRFSource(name, modelType, dict, obr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::fv::solidRotationSource::phiOmega() const
{
    // Zero the flux on any involved boundaries - to compensate for boundary
    // faces that are not exactly normal to the rotational flux
    tmp<surfaceScalarField> po(-GRFSource::phiOmega());
    surfaceScalarField::Boundary& pb = po->boundaryFieldRef();
    forAll(pb, patchi)
    {
        if (!pb[patchi].coupled())
        {
            pb[patchi] *= (1-this->fvmMask_.boundaryField()[patchi]);
        }
    }
    return po;
}

void Foam::fv::solidRotationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    GRFSource::addSup(rho, eqn, fieldI);

    surfaceScalarField rhof(fvc::interpolate(rho));
    volScalarField divPhi(-fvc::div(-1*rhof*phiOmega()));
    // Using SuSp (as done by bounded scheme) leads to slow convergence in
    // diffusion-dominated system as it does not keep source and sink matched
    eqn -= fvm::Sp(divPhi, eqn.psi());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
