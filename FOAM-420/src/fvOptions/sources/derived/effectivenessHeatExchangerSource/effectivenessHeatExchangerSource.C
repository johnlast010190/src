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
    (c) 2013-2015 OpenFOAM Foundation
    (c) 2010-2017, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "effectivenessHeatExchangerSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"
#include "meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "primitives/strings/keyType/keyType.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(effectivenessHeatExchangerSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        effectivenessHeatExchangerSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchangerSource::init()
{
    const faceZone& fZone = mesh_.faceZones()[zoneID_];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        label faceI = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(faceI))
        {
            faceId = faceI;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(faceI);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = faceI - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }
    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);

    calculateTotalArea(faceZoneArea_);

    //find the upwind cells for primaryInletT
    determineUpwindCells();
}


void Foam::fv::effectivenessHeatExchangerSource::calculateTotalArea
(
    scalar& area
)
{
    area = 0;
    forAll(faceId_, i)
    {
        label faceI = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchI = facePatchId_[i];
            area += mesh_.magSf().boundaryField()[patchI][faceI];
        }
        else
        {
            area += mesh_.magSf()[faceI];
        }
    }
    reduce(area, sumOp<scalar>());
}

Foam::labelList Foam::fv::effectivenessHeatExchangerSource::zoneBoundaryFaces() const
{

    DynamicList<label> zoneBoundFaces(0);

    boolList zoneCells(mesh_.nCells(), false);

    forAll(cells_, i)
    {
        zoneCells[cells_[i]] = true;
    }

    const labelList& nei = mesh_.neighbour();
    const labelList& own = mesh_.faceOwner();

    forAll(nei, i)
    {
        if (zoneCells[own[i]] && !zoneCells[nei[i]])
        {
            zoneBoundFaces.append(i);
        }
        else if (zoneCells[nei[i]] && !zoneCells[own[i]])
        {
            zoneBoundFaces.append(i);
        }
    }

    // check on processor patches
    boolList neiZoneID(mesh_.nFaces()-mesh_.nInternalFaces(),false);

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        if (zoneCells[own[faceI]])
        {
            neiZoneID[faceI-mesh_.nInternalFaces()] =  true;
        }
    }

    boolList swappedZoneID(neiZoneID);

    syncTools::swapBoundaryFaceList(mesh_, swappedZoneID);

    const polyPatchList& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
        const processorPolyPatch& procPatch
            = refCast<const processorPolyPatch>(patches[patchI]);

            forAll(procPatch, i)
            {
                label meshFaceI = procPatch.start() + i;
                label bFaceI =  meshFaceI - mesh_.nInternalFaces();

                if
                (
                    !swappedZoneID[bFaceI] && neiZoneID[bFaceI]
                )
                {
                    zoneBoundFaces.append(meshFaceI);
                }
            }
        }
        else
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            forAll(patch, i)
            {
                label meshFaceI = patch.start() + i;
                label bFaceI =  meshFaceI - mesh_.nInternalFaces();
                if (neiZoneID[bFaceI])
                {
                    zoneBoundFaces.append(meshFaceI);
                }
            }
        }
    }

    return labelList(zoneBoundFaces, true);
}

void Foam::fv::effectivenessHeatExchangerSource::determineUpwindCells()
{
    const labelList& nei = mesh_.neighbour();
    const labelList& own = mesh_.faceOwner();

    boolList zoneCells(mesh_.nCells(), false);

    forAll(cells_, i)
    {
        zoneCells[cells_[i]] = true;
    }

    upwindCells_.setSize(faceId_.size(), -1);

    forAll(faceId_, i)
    {
        upwindCells_[i] = (zoneCells[own[faceId_[i]]]) ? nei[faceId_[i]] : own[faceId_[i]];
    }

    // find upwind cell's total volume
    const scalarField& V = mesh_.V();

    forAll(upwindCells_,j)
    {
        inletCellVolume_ += V[upwindCells_[j]];
    }

    reduce(inletCellVolume_, sumOp<scalar>());
}

Foam::tmp<Foam::volScalarField>
Foam::fv::effectivenessHeatExchangerSource::getRho() const
{
    const objectRegistry& obr_ = mesh_;

    if (obr_.found("rho"))
    {
        return(obr_.lookupObject<volScalarField>("rho"));
    }
    else
    {
        dimensionedScalar rho
        (
            "rho",
            obr_.lookupObject<dictionary>("transportProperties")
                .lookup("rho")
        );

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    obr_
                ),
                mesh_,
                rho
            )
        );
    }
}

Foam::scalar Foam::fv::effectivenessHeatExchangerSource::zoneMassFlux
(
    const fvMesh& mesh,
    const labelList& faces
)
{
    scalar fsum = 0;
    const surfaceScalarField& phi
         = mesh.lookupObject<surfaceScalarField>("phi");

    surfaceScalarField rhof(fvc::interpolate(getRho()));

    bool incomp = false;
    if (phi.dimensions() == dimVelocity * dimArea)
    {
        incomp = true;
    }

    forAll(faces, i)
    {
       label meshFaceI = faces[i];
       if (meshFaceI < mesh.nInternalFaces())
       {
           if (incomp)
           {
               fsum += rhof[meshFaceI]*mag(phi[meshFaceI]);
           }
           else
           {
               fsum += mag(phi[meshFaceI]);
           }
       }
       else
       {
          label patchI = mesh.boundaryMesh().whichPatch(meshFaceI);
          const polyPatch& patch = mesh.boundaryMesh()[patchI];
          if (!isA<emptyPolyPatch>(patch))
          {
              label bfI = meshFaceI - patch.start();

              if (incomp)
              {
                  fsum += rhof.boundaryField()[patchI][bfI]
                      *mag(phi.boundaryField()[patchI][bfI]);
              }
              else
              {
                   fsum += mag(phi.boundaryField()[patchI][bfI]);
              }
          }
       }
    }
    reduce(fsum, sumOp<scalar>());
    return fsum;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::effectivenessHeatExchangerSource::effectivenessHeatExchangerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    secondaryMassFlowRate_(readScalar(coeffs_.lookup("secondaryMassFlowRate"))),
    secondaryInletT_(readScalar(coeffs_.lookup("secondaryInletT"))),
    secondaryCp_(readScalar(coeffs_.lookup("secondaryCp"))),
    primaryInletT_(SMALL),
    eTable_(),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    TName_(coeffs_.lookupOrDefault<word>("TName", "T")),
    phiName_(coeffs_.lookupOrDefault<word>("phiName", "phi")),
    faceZoneName_(coeffs_.lookup("faceZone")),
    zoneID_(mesh_.faceZones().findZoneID(faceZoneName_)),
    faceId_(),
    facePatchId_(),
    faceSign_(),
    faceZoneArea_(0),
    upwindCells_(),
    inletCellVolume_(SMALL),
    zoneBoundaryFaces_()
{
    if (zoneID_ < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::effectivenessHeatExchangerSource::initialise()
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const surfaceScalarField& phi
    = mesh_.lookupObject<surfaceScalarField>("phi");

    bool incomp = false;
    if (phi.dimensions() == dimVelocity * dimArea)
    {
        incomp = true;
    }

    if (!incomp)
    {
        const basicThermo& thermo =
            obr_.lookupObject<basicThermo>(basicThermo::dictName);
        fieldNames_.setSize(1, thermo.he().name());
    }
    else
    {
        fieldNames_.setSize(1,TName_);
    }

    applied_.setSize(1, false);

    eTable_.reset(new interpolation2DTable<scalar>(coeffs_));

    init();

    zoneBoundaryFaces_ = zoneBoundaryFaces();

    return true;
}


void Foam::fv::effectivenessHeatExchangerSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

    // getting Cp from the turbulence model
    const incompressible::turbulenceModel & turbModel =
    obr_.lookupObject<incompressible::turbulenceModel>
    (
        turbulenceModel::propertiesName
        );

    const surfaceScalarField Cpf(fvc::interpolate(turbModel.transport().Cp()));

    scalar totalPhi = 0;
    scalar CpfMean = 0;
    forAll(faceId_, i)
    {
        label faceI = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchI = facePatchId_[i];

            CpfMean +=
                Cpf.boundaryField()[patchI][faceI]
        *mesh_.magSf().boundaryField()[patchI][faceI];
        }
        else
        {
            CpfMean += Cpf[faceI]*mesh_.magSf()[faceI];
        }
    }

    reduce(CpfMean, sumOp<scalar>());

    // find primary inlet average T instead of user defined

    primaryInletT_ = 0;

    forAll(upwindCells_, cI)
    {
    primaryInletT_ += T[upwindCells_[cI]]*mesh_.V()[upwindCells_[cI]];
    }

    reduce(primaryInletT_, sumOp<scalar>());
    primaryInletT_ = primaryInletT_/inletCellVolume_;

    //more accurately find mass flow rate
    totalPhi = 0.5 * this->zoneMassFlux(mesh_, zoneBoundaryFaces_);

    scalar Qt =
    eTable_()(mag(totalPhi), secondaryMassFlowRate_)
    *(secondaryInletT_ - primaryInletT_)
    *(CpfMean/faceZoneArea_)*mag(totalPhi);

    const scalarField TCells(T, cells_);
    scalar Tref = 0;
    if (Qt > 0)
    {
        Tref = max(TCells);
        reduce(Tref, maxOp<scalar>());
    }
    else
    {
        Tref = min(TCells);
        reduce(Tref, minOp<scalar>());
    }

    scalarField deltaTCells(cells_.size(), 0);
    forAll(deltaTCells, i)
    {
        if (Qt > 0)
        {
            deltaTCells[i] = max(Tref - TCells[i], 0.0);
        }
        else
        {
            deltaTCells[i] = max(TCells[i] - Tref, 0.0);
        }
    }

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();
    const volScalarField rho( getRho() );
    scalar sumWeight = 0;
    forAll(cells_, i)
    {
        sumWeight += V[cells_[i]]*mag(U[cells_[i]])*deltaTCells[i];
    }

    reduce(sumWeight, sumOp<scalar>());
    sumWeight = max(sumWeight,VSMALL);

    if ((this->V() > VSMALL) && (mag(Qt) > VSMALL) && (sumWeight > VSMALL))
    {
        scalarField& heSource = eqn.source();
        forAll(cells_, i)
        {
        //Teqn for ico solvers -> he/(rho*Cp)
            heSource[cells_[i]] -=
                (Qt/((CpfMean/faceZoneArea_)*rho[cells_[i]]))
        *V[cells_[i]]*mag(U[cells_[i]])*deltaTCells[i]/sumWeight;
        }
    }

    scalar secondaryOutletT(0);


    if (secondaryMassFlowRate_ > VSMALL)
    {
    secondaryOutletT = secondaryInletT_ - Qt/( secondaryCp_*secondaryMassFlowRate_);

    }

    if (Pstream::master())
    {
        Info<< indent << "\nfvOption:effectivenessHeatExchangerSource:" << this->name() << nl;
        Info<< indent << "  Net mass flux [Kg/s] = " << totalPhi << nl;
        Info<< indent << "  Total energy exchange [W] = " << Qt << nl;
        Info<< indent << "  TPrimaryInlet [K] = " << primaryInletT_ << nl;
    Info<< indent << "  TSecondaryOutlet [K] = " << secondaryOutletT << nl;
        Info<< indent << "  Tref [K] = " << Tref << nl;
        Info<< indent << "  Efficiency : "
            << eTable_()(mag(totalPhi), secondaryMassFlowRate_) << nl << nl;
    }
}

void Foam::fv::effectivenessHeatExchangerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label
)
{
    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);
    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

    const surfaceScalarField Cpf(fvc::interpolate(thermo.Cp()));

    scalar totalPhi = 0;
    scalar CpfMean = 0;
    forAll(faceId_, i)
    {
        label faceI = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchI = facePatchId_[i];

            CpfMean +=
                Cpf.boundaryField()[patchI][faceI]
               *mesh_.magSf().boundaryField()[patchI][faceI];
        }
        else
        {
            CpfMean += Cpf[faceI]*mesh_.magSf()[faceI];
        }
    }

    reduce(CpfMean, sumOp<scalar>());

    // find primary inlet average T instead of user defined

    primaryInletT_ = 0;

    forAll(upwindCells_, cI)
    {
        primaryInletT_ += T[upwindCells_[cI]]*mesh_.V()[upwindCells_[cI]];
    }

    reduce(primaryInletT_, sumOp<scalar>());
    primaryInletT_ = primaryInletT_/inletCellVolume_;

    //more accurately find mass flow rate
    totalPhi = 0.5 * this->zoneMassFlux(mesh_, zoneBoundaryFaces_);

    scalar Qt =
        eTable_()(mag(totalPhi), secondaryMassFlowRate_)
       *(secondaryInletT_ - primaryInletT_)
       *(CpfMean/faceZoneArea_)*mag(totalPhi);

    const scalarField TCells(T, cells_);
    scalar Tref = 0;
    if (Qt > 0)
    {
        Tref = max(TCells);
        reduce(Tref, maxOp<scalar>());
    }
    else
    {
        Tref = min(TCells);
        reduce(Tref, minOp<scalar>());
    }

    scalarField deltaTCells(cells_.size(), 0);
    forAll(deltaTCells, i)
    {
        if (Qt > 0)
        {
            deltaTCells[i] = max(Tref - TCells[i], 0.0);
        }
        else
        {
            deltaTCells[i] = max(TCells[i] - Tref, 0.0);
        }
    }

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();
    scalar sumWeight = 0;
    forAll(cells_, i)
    {
        sumWeight += V[cells_[i]]*mag(U[cells_[i]])*deltaTCells[i];
    }

    reduce(sumWeight, sumOp<scalar>());
    sumWeight = max(sumWeight,VSMALL);

    if ((this->V() > VSMALL) && (mag(Qt) > VSMALL) && (sumWeight > VSMALL))
    {
        scalarField& heSource = eqn.source();

        forAll(cells_, i)
        {
            heSource[cells_[i]] -=
                Qt*V[cells_[i]]*mag(U[cells_[i]])*deltaTCells[i]/sumWeight;
        }
    }

    scalar secondaryOutletT(0);

    if (secondaryMassFlowRate_ > VSMALL)
    {
        secondaryOutletT = secondaryInletT_ - Qt/( secondaryCp_*secondaryMassFlowRate_);

    }


    if (Pstream::master())
    {
        Info<< indent << "\nfvOption:effectivenessHeatExchangerSource:" << this->name() << nl;
        Info<< indent << "  Net mass flux [Kg/s] = " << totalPhi << nl;
        Info<< indent << "  Total energy exchange [W] = " << Qt << nl;
        Info<< indent << "  TPrimaryInlet [K] = " << primaryInletT_ << nl;
        Info<< indent << "  TSecondaryOutlet [K] = " << secondaryOutletT << nl;
        Info<< indent << "  Tref [K] = " << Tref << nl;
        Info<< indent << "  Efficiency : "
            << eTable_()(mag(totalPhi), secondaryMassFlowRate_) << nl << nl;
    }
}


bool Foam::fv::effectivenessHeatExchangerSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("secondaryMassFlowRate") >> secondaryMassFlowRate_;
        coeffs_.lookup("secondaryInletT") >> secondaryInletT_;
        coeffs_.lookup("secondaryCp") >> secondaryCp_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
