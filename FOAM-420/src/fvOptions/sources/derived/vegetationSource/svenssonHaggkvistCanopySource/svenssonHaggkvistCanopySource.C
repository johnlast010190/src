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

#include "svenssonHaggkvistCanopySource.H"
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
    defineTypeNameAndDebug(svenssonHaggkvistCanopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        svenssonHaggkvistCanopySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::svenssonHaggkvistCanopySource::svenssonHaggkvistCanopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
  vegetationSource(name, modelType, dict, obr),
  CpEps1_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "CpEps1",
            coeffs_,
            1.8
        )
  )
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::fv::svenssonHaggkvistCanopySource::addSup
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


void Foam::fv::svenssonHaggkvistCanopySource::addSup
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


void Foam::fv::svenssonHaggkvistCanopySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
  const volScalarField& canopy = canopy_;

  if (eqn.psi().name() == word("k"))
  {
      const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
      eqn += canopy*pow(mag(U),3);
  }
  else if (eqn.psi().name() == word("epsilon"))
  {
      const volScalarField& epsilon = eqn.psi();
      const volScalarField& k = obr_.lookupObject<volScalarField>("k");
      const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

      fvMatrix<scalar> Sepsilon
      (
      fvm::Sp(canopy*CpEps1_/k*pow(mag(U),3), epsilon)
      );

      eqn += Sepsilon;
  }
  else if (eqn.psi().name() == word("omega"))
  {
      const volScalarField& omega = eqn.psi();
      const volScalarField& k = obr_.lookupObject<volScalarField>("k");
      const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

      fvMatrix<scalar> Somega
      (
      fvm::Sp(canopy*CpEps1_/k*pow(mag(U),3), omega)
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


void Foam::fv::svenssonHaggkvistCanopySource::addSup
 (
     const volScalarField& rho,
     fvMatrix<scalar>& eqn,
     const label fieldi
 )
 {

  const volScalarField& canopy = canopy_;

  if (eqn.psi().name() == word("k"))
  {
      const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

      eqn += rho*canopy*pow(mag(U),3);
  }
  else if (eqn.psi().name() == word("epsilon"))
  {

      const volScalarField& epsilon = eqn.psi();
      const volScalarField& k = obr_.lookupObject<volScalarField>("k");
      const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

      fvMatrix<scalar> Sepsilon
      (
      fvm::Sp(rho*canopy*CpEps1_/k*pow(mag(U),3), epsilon)
      );

    eqn += Sepsilon;
  }
  else if (eqn.psi().name() == word("omega"))
  {

      const volScalarField& omega = eqn.psi();
      const volScalarField& k = obr_.lookupObject<volScalarField>("k");
      const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

      fvMatrix<scalar> Somega
      (
      fvm::Sp(rho*canopy*CpEps1_/k*pow(mag(U),3), omega)
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
