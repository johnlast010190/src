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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "dalpeMassonCanopySource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/fvOptions/fvOption.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(dalpeMassonCanopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        dalpeMassonCanopySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::dalpeMassonCanopySource::dalpeMassonCanopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
  vegetationSource(name, modelType, dict, obr),
  betaP_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaP",
            coeffs_,
            1.0
        )
  ),
  betaD_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaD",
            coeffs_,
            5.03
        )
  ),
  C4_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "C4",
            coeffs_,
            0.78
        )
  ),
  C5_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "C5",
            coeffs_,
            0.78
        )
   )
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::dalpeMassonCanopySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{

    if (eqn.psi().name() == word(UName_))
    {

        const volVectorField& U = eqn.psi();
        const volScalarField& canopy = canopy_;

        fvMatrix<vector> S_canopy
        (
            fvm::Sp(canopy * mag(U), U)
        );

        eqn -=  S_canopy;
    }
}


void Foam::fv::dalpeMassonCanopySource::addSup
 (
     const volScalarField& rho,
     fvMatrix<vector>& eqn,
     const label fieldi
 )
{

    if (eqn.psi().name() == word(UName_))
    {

        const volVectorField& U = eqn.psi();
        const volScalarField& canopy = canopy_;

        fvMatrix<vector> S_canopy
        (
            fvm::Sp(rho*canopy * mag(U), U)
        );

        eqn -=  S_canopy;
    }
}


void Foam::fv::dalpeMassonCanopySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{

    const volScalarField& canopy = canopy_;

    if (eqn.psi().name() == word("k"))
    {

        const volScalarField& k = eqn.psi();
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        fvMatrix<scalar> Sk
        (
            betaP_*canopy*pow(mag(U),3) - fvm::Sp(betaD_*canopy*mag(U), k)
        );

        eqn +=  Sk;
    }
    else if (eqn.psi().name() == word("epsilon"))
    {

        const volScalarField& epsilon = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        fvMatrix<scalar> Sepsilon
        (
            fvm::Sp(canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), epsilon)
        );

        eqn += Sepsilon;

    }
    else if (eqn.psi().name() == word("omega"))
    {

        const volScalarField& omega = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        // Somega = Sepsilon/k
        fvMatrix<scalar> Somega
        (
            fvm::Sp(canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), omega)
        );

        eqn += Somega;

    }
    else if (eqn.psi().name() == word(TName_))
    {
        const scalarField& V = mesh_.V();
        const volScalarField rho( getRho() );
        const volScalarField Cp( getCp() );

        scalarField& TSource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            TSource[celli] -=LAD[i]*Qs[i]/(rho[celli]*Cp[celli])*V[celli];
        }
    }
    else if (eqn.psi().name() == word("w"))
    {
        const scalarField& V = mesh_.V();
        const volScalarField rho( getRho() );

        scalarField& humiditySource( eqn.source() );

        forAll(cells_, i)
        {
            label celli = cells_[i];
            humiditySource[celli] -=LAD[i]*Ql[i]/(rho[celli]*lambdaVap_)*V[celli];
        }
    }

}


void Foam::fv::dalpeMassonCanopySource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{

    const volScalarField& canopy = canopy_;

    if (eqn.psi().name() == word("k"))
    {

        const volScalarField& k = eqn.psi();
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        fvMatrix<scalar> Sk
        (
            betaP_*rho*canopy*pow(mag(U),3) - fvm::Sp(betaD_*rho*canopy*mag(U), k)
        );

        eqn +=  Sk;
    }
    else if (eqn.psi().name() == word("epsilon"))
    {
        const volScalarField& epsilon = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        fvMatrix<scalar> Sepsilon
        (
            fvm::Sp(rho*canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), epsilon)
        );

        eqn += Sepsilon;
    }
    else if (eqn.psi().name() == word("omega"))
    {
        const volScalarField& omega = eqn.psi();
        const volScalarField& k = obr_.lookupObject<volScalarField>("k");
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        // Somega = Sepsilon/k
        fvMatrix<scalar> Somega
        (
            fvm::Sp(rho*canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), omega)
        );

        eqn += Somega;
    }
    else if (eqn.psi().name() == word(eName_))
    {
        const scalarField& V = mesh_.V();

        scalarField& heSource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            heSource[celli] -=LAD[i]*Qs[i]*V[celli];
        }
    }
    else if (eqn.psi().name() == word("w"))
    {
        const scalarField& V = mesh_.V();

        scalarField& humiditySource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            humiditySource[celli] -=LAD[i]*Ql[i]/lambdaVap_*V[celli];
        }
    }

}



// ************************************************************************* //
