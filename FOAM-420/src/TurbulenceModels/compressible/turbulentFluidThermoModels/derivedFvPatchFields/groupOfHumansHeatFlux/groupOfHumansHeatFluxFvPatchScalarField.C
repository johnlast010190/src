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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "groupOfHumansHeatFluxFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "patchDist/wallPointData/wallPointData.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<Foam::compressible::groupOfHumansHeatFluxFvPatchScalarField::eGroup, 4>::names[] =
{
    "maleAdult",
    "femaleAdult",
    "maleChild",
    "femaleChild"
};

namespace Foam
{
namespace compressible
{
// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


const Foam::NamedEnum
    <Foam::compressible::groupOfHumansHeatFluxFvPatchScalarField::eGroup, 4>
    Foam::compressible::groupOfHumansHeatFluxFvPatchScalarField::eGroupNames_;


// * * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * //

void groupOfHumansHeatFluxFvPatchScalarField::findLocalCells
(
    const fvMesh& mesh,
    scalar d
)
{
    // find cells with distance from pacth smaller than d

    label nBoundaryFaces = this->size();
    DynamicList<wallPointData<label>> faceDist(nBoundaryFaces);
    DynamicList<label> changedFaces(nBoundaryFaces);

    forAll(patch().Cf(), patchFaceI)
    {
        label meshFaceI = patch().patch().start() + patchFaceI;

        changedFaces.append(meshFaceI);

        faceDist.append
        (
            wallPointData<label>
            (
                patch().Cf()[patchFaceI],
                patch().index(), // passive label
                0.0
            )
        );
    }

    changedFaces.shrink();
    faceDist.shrink();

    // Perform mesh wave to calculate distance to patch
    List<wallPointData<label>> faceInfo(mesh.nFaces()), cellInfo(mesh.nCells());
    FaceCellWave<wallPointData<label>> waveInfo
    (
        mesh,
        changedFaces,
        faceDist,
        faceInfo,
        cellInfo,
        mesh.globalData().nTotalCells()    // max iterations
    );

    // Cells to cellset nearCells_
    scalar maxDsqr = sqr(d);
    nearVolume_ = 0;
    forAll(mesh.cells(), cellI)
    {
        label patchI = cellInfo[cellI].data();
        scalar dSqr = cellInfo[cellI].distSqr();

        if (patchI != -1 && dSqr < maxDsqr)
        {
            nearCells_.insert(cellI);
            nearVolume_ += mesh.V()[cellI];
        }
    }

    reduce(nearVolume_, sumOp<scalar>());
    Areal_ = gSum(patch().magSf());
}

scalar groupOfHumansHeatFluxFvPatchScalarField::humanSurfaceArea
(
    scalar W,
    scalar H
) const
{
    return (0.202 * ::pow(W, 0.425) * ::pow(H, 0.725));
}

void groupOfHumansHeatFluxFvPatchScalarField::staticStatisticsCalculator()
{
    // calculate effective values
    Ntotal_ = 0;
    Weff_ = 0;
    Heff_ = 0;

    forAll(people_, gI)
    {
        Weff_ += people_[gI].N * people_[gI].weight;
        Heff_ += people_[gI].N * people_[gI].height;
        Ntotal_ += people_[gI].N;
    }

    Weff_ /= Ntotal_;
    Heff_ /= Ntotal_;

    Aeff_ = humanSurfaceArea(Weff_, Heff_);
}

scalar groupOfHumansHeatFluxFvPatchScalarField::humanHeatFlux() const
{
    scalar adultMaleSurfaceArea = humanSurfaceArea
    (
        people_[mA].weight, people_[mA].height
    );

    // scale the heat flux with the supposed area to get the right overall
    // heat release rate
    scalar Atotal = 0;
    forAll(people_, gI)
    {
        Atotal
            += people_[gI].N
            * humanSurfaceArea(people_[gI].weight, people_[gI].height);
    }

    return
    (
        (-0.875 * (Troom() - 273.15) + 138.5)
        *Atotal / Areal_ / adultMaleSurfaceArea
    );
}

scalar groupOfHumansHeatFluxFvPatchScalarField::effectiveMetabolicRate() const
{
    return
    (
        Mstd_ * (people_[mA].N + 0.9*people_[fA].N
        + (-1.39*(people_[mC].age - 5) + 60.0)/46.5 * people_[mC].N
        + (-0.88 * (people_[fC].age - 5) + 61.6)/46.5 * people_[fC].N)
        / Ntotal_
    );
}

scalar groupOfHumansHeatFluxFvPatchScalarField::meanRadiantTemperature() const
{
    scalar Tmr = 0;

    if (db().foundObject<volScalarField>("qin"))
    {
        volScalarField Qin(db().lookupObject<volScalarField>("qin"));

        scalarField Tr
        (
            Foam::pow((mag(Qin.boundaryField()[patch().index()])
            / (emissivityFarWall_ * emissivityLocal_ * sigmaSB_))(), 0.25)
            * patch().magSf()
        );

        Tmr = gSum(Tr)/Areal_;
    }

    return Tmr;
}

List<groupOfHumansHeatFluxFvPatchScalarField::demographic>
groupOfHumansHeatFluxFvPatchScalarField::makePeople
(
    const dictionary& dict
) const
{
    List<demographic> people(4);

    forAllConstIter(dictionary, dict, iter)
    {
        const word& nk = iter().keyword();
        const dictionary& dict = iter().dict();
        people[eGroupNames_[nk]] = demographic(dict);
    }

    return people;
}

scalar groupOfHumansHeatFluxFvPatchScalarField::effectiveRoomtemperature() const
{
    scalar Ulocal = calcLocalMean<volVectorField>("U");
    scalar Urel = Ulocal + 0.0052*(Meff() - 58.0);
    scalar Troom = Tlocal();
    scalar TMeanRadiant = meanRadiantTemperature();

    if (TMeanRadiant > 0)
    {
        scalar k = 0;
        if (Urel < 0.2)
        {
            k = 0.5;
        }
        else if (Urel < 0.6)
        {
            k = 0.6;
        }
        else
        {
            k = 0.7;
        }

        Troom = k * Tlocal() + (1 - k ) * TMeanRadiant;
    }

    return Troom;
}

template<class volField>
scalar groupOfHumansHeatFluxFvPatchScalarField::calcLocalMean(word f) const
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    const volField& F = this->db().lookupObject<volField>(f);

    scalar Flocal = 0;

    forAllConstIter(cellSet, nearCells_, iter)
    {
        label cellI = iter.key();

        Flocal += mesh.V()[cellI] * mag(F[cellI]);
    }

    reduce(Flocal, sumOp<scalar>());
    Flocal /= nearVolume_;

    return Flocal;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

groupOfHumansHeatFluxFvPatchScalarField::
groupOfHumansHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(p, iF),
    people_(4),
    Mstd_(0),
    dlocal_(0),
    emissivityFarWall_(0.9),
    emissivityLocal_(0.8),
    sigmaSB_(5.6704e-8),
    Weff_(0),
    Heff_(0),
    Aeff_(0),
    Areal_(0),
    nearCells_(p.patch().boundaryMesh().mesh(), "nearCells", 0),
    nearVolume_(0),
    Ntotal_(0),
    MeffPtr_(),
    PlocalPtr_(),
    TlocalPtr_(),
    TroomPtr_()
{}


groupOfHumansHeatFluxFvPatchScalarField::
groupOfHumansHeatFluxFvPatchScalarField
(
    const groupOfHumansHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    people_(ptf.people_),
    Mstd_(ptf.Mstd_),
    dlocal_(ptf.dlocal_),
    emissivityFarWall_(ptf.emissivityFarWall_),
    emissivityLocal_(ptf.emissivityLocal_),
    sigmaSB_(ptf.sigmaSB_),
    Weff_(ptf.Weff_),
    Heff_(ptf.Heff_),
    Aeff_(ptf.Aeff_),
    Areal_(ptf.Areal_),
    nearCells_
    (
        p.patch().boundaryMesh().mesh(),
        "nearCells",
        p.patch().boundaryMesh().mesh().nCells()
    ),
    nearVolume_(ptf.nearVolume_),
    Ntotal_(ptf.Ntotal_),
    MeffPtr_(),
    PlocalPtr_(),
    TlocalPtr_(),
    TroomPtr_()
{}


groupOfHumansHeatFluxFvPatchScalarField::
groupOfHumansHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(p, iF, dict),
    people_(makePeople(dict.subDict("demographics"))),
    Mstd_(gAverage(turbulentHeatFluxTemperatureFvPatchScalarField::q_)),
    dlocal_(readScalar(dict.lookup("averagingDistance"))),
    emissivityFarWall_(dict.lookupOrDefault<scalar>("emissivityWalls", 0.9)),
    emissivityLocal_(dict.lookupOrDefault<scalar>("emissivityPeople", 0.8)),
    sigmaSB_(dict.lookupOrDefault<scalar>("sigmaSB", 5.6704e-8)),
    Weff_(0),
    Heff_(0),
    Aeff_(0),
    Areal_(0),
    nearCells_(p.patch().boundaryMesh().mesh(), "nearCells", 0),
    nearVolume_(0),
    Ntotal_(0),
    MeffPtr_(),
    PlocalPtr_(),
    TlocalPtr_(),
    TroomPtr_()
{
    //recalculate derived data
    bool solving = true;

    if (!Pstream::parRun())
    {
        const fvMesh& mesh(patch().boundaryMesh().mesh());

        forAll(mesh.boundary(), pI)
        {
            if (mesh.boundary()[pI].coupled())
            {
                solving = false;
                break;
            }
        }
    }


    if (solving)
    {
        findLocalCells(p.boundaryMesh().mesh(), dlocal_);
        staticStatisticsCalculator();
    }
}


groupOfHumansHeatFluxFvPatchScalarField::
groupOfHumansHeatFluxFvPatchScalarField
(
    const groupOfHumansHeatFluxFvPatchScalarField& ptf
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(ptf),
    people_(ptf.people_),
    Mstd_(ptf.Mstd_),
    dlocal_(ptf.dlocal_),
    emissivityFarWall_(ptf.emissivityFarWall_),
    emissivityLocal_(ptf.emissivityLocal_),
    sigmaSB_(ptf.sigmaSB_),
    Weff_(ptf.Weff_),
    Heff_(ptf.Heff_),
    Aeff_(ptf.Aeff_),
    Areal_(ptf.Areal_),
    nearCells_
    (
        ptf.patch().patch().boundaryMesh().mesh(),
        "nearCells",
        ptf.patch().patch().boundaryMesh().mesh().nCells()
    ),
    nearVolume_(ptf.nearVolume_),
    Ntotal_(ptf.Ntotal_),
    MeffPtr_(),
    PlocalPtr_(),
    TlocalPtr_(),
    TroomPtr_()
{

    bool solving = true;

    if (!Pstream::parRun())
    {
        const fvMesh& mesh(patch().boundaryMesh().mesh());

        forAll(mesh.boundary(), pI)
        {
            if (mesh.boundary()[pI].coupled())
            {
                solving = false;
                break;
            }
        }
    }

    if (solving)
    {
        //recalculate derived data
        findLocalCells(refCast<const fvMesh>(this->db()), dlocal_);
        staticStatisticsCalculator();
    }
}


groupOfHumansHeatFluxFvPatchScalarField::
groupOfHumansHeatFluxFvPatchScalarField
(
    const groupOfHumansHeatFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(ptf, iF),
    people_(ptf.people_),
    Mstd_(ptf.Mstd_),
    dlocal_(ptf.dlocal_),
    emissivityFarWall_(ptf.emissivityFarWall_),
    emissivityLocal_(ptf.emissivityLocal_),
    sigmaSB_(ptf.sigmaSB_),
    Weff_(ptf.Weff_),
    Heff_(ptf.Heff_),
    Aeff_(ptf.Aeff_),
    Areal_(ptf.Areal_),
    nearCells_
    (
        ptf.patch().patch().boundaryMesh().mesh(),
        "nearCells",
        ptf.patch().patch().boundaryMesh().mesh().nCells()
    ),
    nearVolume_(ptf.nearVolume_),
    Ntotal_(ptf.Ntotal_),
    MeffPtr_(),
    PlocalPtr_(),
    TlocalPtr_(),
    TroomPtr_()
{

    bool solving = true;

    if (!Pstream::parRun())
    {
        const fvMesh& mesh(patch().boundaryMesh().mesh());

        forAll(mesh.boundary(), pI)
        {
            if (mesh.boundary()[pI].coupled())
            {
                solving = false;
                break;
            }
        }
    }

    if (solving)
    {
        //recalculate derived data
        findLocalCells(ptf.patch().boundaryMesh().mesh(), dlocal_);
        staticStatisticsCalculator();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void groupOfHumansHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    TlocalPtr_.reset(new scalar(calcLocalMean<volScalarField>("T")));

    {
        volScalarField Pabs
        (
            "Pabs",
            this->db().lookupObject<volScalarField>("p")
          + this->db().lookupObject<basicThermo>(basicThermo::dictName).pRef()
        );
        PlocalPtr_.reset(new scalar(calcLocalMean<volScalarField>("Pabs")));
    }

    TroomPtr_.reset(new scalar(effectiveRoomtemperature()));

    Mstd_ = humanHeatFlux();

    turbulentHeatFluxTemperatureFvPatchScalarField::q_ = Mstd_;

    //update Meff
    MeffPtr_.reset(new scalar(effectiveMetabolicRate()));

    turbulentHeatFluxTemperatureFvPatchScalarField::updateCoeffs();
}


scalar groupOfHumansHeatFluxFvPatchScalarField::Plocal() const
{
    if (!PlocalPtr_.valid())
    {
        volScalarField Pabs
        (
            "Pabs",
            this->db().lookupObject<volScalarField>("p")
          + this->db().lookupObject<basicThermo>(basicThermo::dictName).pRef()
        );
        PlocalPtr_.reset(new scalar(calcLocalMean<volScalarField>("Pabs")));
    }

    return PlocalPtr_();
}


scalar groupOfHumansHeatFluxFvPatchScalarField::Tlocal() const
{
    if (!TlocalPtr_.valid())
    {
        TlocalPtr_.reset(new scalar(calcLocalMean<volScalarField>("T")));
    }

    return TlocalPtr_();
}


scalar groupOfHumansHeatFluxFvPatchScalarField::Troom() const
{
    if (!TroomPtr_.valid())
    {
        TroomPtr_.reset(new scalar(effectiveRoomtemperature()));
    }

    return TroomPtr_();
}


scalar groupOfHumansHeatFluxFvPatchScalarField::Meff() const
{
    if (!MeffPtr_.valid())
    {
        MeffPtr_.reset(new scalar(effectiveMetabolicRate()));
    }

    return MeffPtr_();
}


void groupOfHumansHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    dictionary demographics("demographics");
    forAll(people_, gI)
    {
        dictionary dict
        (
            groupOfHumansHeatFluxFvPatchScalarField::eGroupNames_.names[gI]
        );

        people_[gI].addTo(dict);

        demographics.add(word(dict.name()), dict);
    }

    os.writeKeyword("demographics") << demographics << nl;

    os.writeEntry("averagingDistance", dlocal_);

    writeEntryIfDifferent<scalar>
        (os, "emissivityWalls", 0.9, emissivityFarWall_);
    writeEntryIfDifferent<scalar>
        (os, "emissivityPeople", 0.8, emissivityLocal_);
    writeEntryIfDifferent<scalar>(os, "sigmaSB", 5.6704e-8, sigmaSB_);

    turbulentHeatFluxTemperatureFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    groupOfHumansHeatFluxFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
