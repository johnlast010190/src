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

#include "disperseEulerianSource.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(disperseEulerianSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        disperseEulerianSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::disperseEulerianSource::disperseEulerianSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr),
    des_(obr_.lookupObject<disperseEulerian>(word(dict.lookup("dispersePhase")))),
    heatTrans_(),
    massTrans_(),
    forcesD_(des_.dragModel()),
    forcesL_(des_.liftModel()),
    forcesTD_(des_.tdModel()),
    disableBoundaryForces_(coeffs_.lookupOrDefault<bool>("disableBoundaryForces", false)),
    patches_(coeffs_.lookupOrDefault<wordList>("patches", mesh_.boundaryMesh().names())) //wordList(0)
{
    coeffs_.lookup("fields") >> fieldNames_;

    if (disableBoundaryForces_)
    {
        Info<<"disable forces : " << disableBoundaryForces_
            <<" near boundary patches : " << patches_ << endl;
    }

/*
    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }
*/
    IOdictionary disperseEulerianProperties
    (
        IOobject
        (
            "disperseEulerianProperties",
            mesh_.time().constant(),
            obr_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    PtrList<entry> phaseEntries(disperseEulerianProperties.lookup("phases"));
    heatTrans_.setSize(des_.phases().size());
    massTrans_.setSize(des_.phases().size());

    forAll(heatTrans_, phaseI)
    {
        heatTrans_.set
        (
            phaseI,
            decoupledEulerian::heatTransferModel::New
            (
                phaseEntries[phaseI].dict(),
                des_.phases()[phaseI]
            )
        );

        massTrans_.set
        (
            phaseI,
            decoupledEulerian::massTransferModel::New
            (
                phaseEntries[phaseI].dict(),
                des_.phases()[phaseI]
            )
        );
    }

    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::disperseEulerianSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == "T")
    {
        if (debug)
        {
            Info<<"Inserting thermal source terms into transport Equation for field " << eqn.psi().name() << endl;
        }

        const PtrList<decoupledEulerian::phase>& phases = des_.phases();
        const volScalarField& Tc = eqn.psi();

        for (label i=0; i < phases.size(); i++)
        {
            volScalarField Kh( heatTrans_[i].Kh() );
            if (debug)
            {
                Info<< "min/max(Kh): " << min(Kh).value() << " , " << max(Kh).value() << endl;
                Info<< "min/max(Tc): " << min(Tc).value() << " , " << max(Tc).value() << endl;
                Info<< "Td: " << phases[i].Td().value() << endl;
            }
            eqn += Kh * phases[i].Td() - fvm::Sp(Kh, Tc) + massTrans_[i].mDot() * phases[i].dHevap() / (max(1.-phases[i].alphad(), scalar(0.01)) * phases[i].rhoc());
        }
    }
    else if (std::string(eqn.psi().name()).find("alpha") != std::string::npos)
    {
        const PtrList<decoupledEulerian::phase>& phases = des_.phases();
        for (label i=0; i < phases.size(); i++)
        {
            if (phases[i].alphad().name() == eqn.psi().name())
            {
                eqn -= massTrans_[i].mDot()/phases[i].rhod();

                if (debug)
                {
                    Info<<"Inserting mass source term into transport Equation for field " << eqn.psi().name() << endl;
                    Info<<"min mass source: " << min(massTrans_[i].mDot()/phases[i].rhod()) << "\n"
                        <<"max mass source: " << max(massTrans_[i].mDot()/phases[i].rhod()) << endl;
                }
            }
        }
    }
}


void Foam::fv::disperseEulerianSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == "U")
    {
        if (debug)
        {
            Info<<"Inserting momentum source terms into transport Equation for field " << eqn.psi().name() << endl;
        }

        const PtrList<decoupledEulerian::phase>& phases = des_.phases();
        const volVectorField& U = eqn.psi();

        for (label i=0; i < phases.size(); i++)
        {
            const volVectorField& Ud = phases[i].Ud();
            const volScalarField scaleFactor( phases[i].alphad() * phases[i].rhod() / phases[i].rhoc() );
            volScalarField Kd(forcesD_[i].K() * scaleFactor);

            volVectorField liftForces( forcesL_[i].F()*scaleFactor );
            volVectorField turbulentDispersionForces( forcesTD_[i].F()*scaleFactor );

            if (disableBoundaryForces_)
            {
                forAll(patches_, pI)
                {
                    word substring = patches_[pI];

                    forAll(mesh_.boundaryMesh(), patchI)
                    {
                        const polyPatch& patch(mesh_.boundaryMesh()[patchI]);
                        word patchName(patch.name());

                        if (patchName.find(substring, 0) != string::npos)
                        {
                            Info<<"found patch " << patchName << " in list " << patches_ << endl;
                            const labelList& cells = mesh_.boundaryMesh()[patchI].faceCells();
                            Info<< "list of cells : " << cells.size() << endl;
                            forAll(cells, cI)
                            {
                                Kd[cI] = 0.0;
                                liftForces[cI] = vector::zero;
                                turbulentDispersionForces[cI] = vector::zero;
                            }
                        }
                    }
                }
            }

            eqn -=
            (
                fvm::Sp(Kd, U)
              - Kd * Ud
              //  Kd * (U - Ud)
              + liftForces
              - turbulentDispersionForces
            );
        }
    }
}


void Foam::fv::disperseEulerianSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == "U")
    {
        if (debug)
        {
            Info<<"Inserting momentum source terms into transport Equation for field " << eqn.psi().name() << endl;
        }

        const PtrList<decoupledEulerian::phase>& phases = des_.phases();
        const volVectorField& U = eqn.psi();

        for (label i=0; i < phases.size(); i++)
        {
            const volVectorField& Ud = phases[i].Ud();
            const volScalarField scaleFactor( phases[i].alphad() * phases[i].rhod() );
            volScalarField Kd(forcesD_[i].K() * scaleFactor);

            volVectorField liftForces( forcesL_[i].F()*scaleFactor );
            volVectorField turbulentDispersionForces( forcesTD_[i].F()*scaleFactor );

            if (disableBoundaryForces_)
            {
                forAll(patches_, pI)
                {
                    word substring = patches_[pI];

                    forAll(mesh_.boundaryMesh(), patchI)
                    {
                        const polyPatch& patch(mesh_.boundaryMesh()[patchI]);
                        word patchName(patch.name());

                        if (patchName.find(substring, 0) != string::npos)
                        {
                            Info<<"found patch " << patchName << " in list " << patches_ << endl;
                            const labelList& cells = mesh_.boundaryMesh()[patchI].faceCells();
                            Info<< "list of cells : " << cells.size() << endl;
                            forAll(cells, cI)
                            {
                                Kd[cI] = 0.0;
                                liftForces[cI] = vector::zero;
                                turbulentDispersionForces[cI] = vector::zero;
                            }
                        }
                    }
                }
            }

            eqn -=
            (
                fvm::Sp(Kd, U)
              - Kd * Ud
              //  Kd * (U - Ud)
              + liftForces
              - turbulentDispersionForces
            );
        }
    }
}


// ************************************************************************* //
