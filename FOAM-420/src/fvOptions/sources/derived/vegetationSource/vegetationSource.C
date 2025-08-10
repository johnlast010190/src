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
    (c) 2018-2019 Esi Ltd.
    (c) 2017 D.Segersson SMHI
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "vegetationSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "transportModel/transportModel.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(vegetationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        vegetationSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::vegetationSource::vegetationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    origin_(coeffs_.lookupOrDefault<vector>("origin", vector::zero)),
    zoneDir_(coeffs_.lookupOrDefault<vector>("direction",vector(1, 0, 0))),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    TName_(coeffs_.lookupOrDefault<word>("TName", "T")),
    PName_(coeffs_.lookupOrDefault<word>("PName", "p")),
    eName_(coeffs_.lookupOrDefault<word>("eName", "h")),
    uniformLAD_(coeffs_.lookupOrDefault<bool>("uniformLAD", false)),
    height_(0),
    Cd_(coeffs_.lookupOrDefault<scalar>("Cd", 0.2)),
    lsize_(coeffs_.lookupOrDefault<scalar>("leafSize", 0.1)),
    LAI_(coeffs_.lookupOrDefault<scalar>("LAI", 0.0)),
    LADmax_(coeffs_.lookupOrDefault<scalar>("LADmax", -1.0)),
    as_(coeffs_.lookupOrDefault<scalar>("as", 0.5)),
    at_(coeffs_.lookupOrDefault<scalar>("at", 0.04)),
    Mvap_(18.02),
    Mair_(28.96),
    lambdaVap_(2500000),
    LADProfile_(coeffs_.lookupOrDefault<scalarList>("LADProfile", scalarList(1, 1.) )),
    LAD(cells_.size(), 0),
    Qs(cells_.size(), 0),
    Ql(cells_.size(), 0)
{
    if ((!uniformLAD_) && (LADProfile_.size() == 1))
    {
        FatalErrorInFunction
            << "Non uniform LAD. LADProfile must be specified."
            << abort(FatalError);
    }

    coeffs_.lookup("fields") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating vegetation source zone: "
        << this->name() << endl;

    if (mag(zoneDir_) < SMALL)
    {
        FatalErrorInFunction
            << "Badly defined axis: zero magnitude: " << zoneDir_
            << abort(FatalError);
    }
    // normalise axis
    zoneDir_ /= mag(zoneDir_);

    //Calculate cellZone height
    PackedBoolList sourcePts(mesh_.nPoints(), 0);

    forAll(cells_, i)
    {
        label celli = cells_[i];

        const labelList& cellPoints = mesh_.cellPoints()[celli];
        forAll(cellPoints, cPtI)
        {
            sourcePts.set(cellPoints[cPtI], 1);
        }
    }

    scalar minHeight = GREAT;
    scalar maxHeight = -GREAT;
    forAll(mesh_.points(), pointi)
    {
        if (sourcePts.get(pointi) == 1)
        {
            vector v(mesh_.points()[pointi] - origin_);
            scalar parallel = (v & zoneDir_);
            minHeight = min(minHeight, parallel);
            maxHeight = max(maxHeight, parallel);
        }
    }
    reduce(minHeight, minOp<scalar>());
    reduce(maxHeight, maxOp<scalar>());

    height_ = mag(maxHeight - minHeight);
    scalar zoneSection_=V_/height_;

    Info<<"    - Canopy height [m]: "<< height_ <<endl;
    Info<<"    - Canopy cross section [m2]: "<< zoneSection_ <<endl;

    init();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::vegetationSource::init()
{

    // if LADmax_ not specified, this will be calculated from LAI
    if (LADmax_ == -1.0 && LAI_ != 0.0)
    {
        Info<<" LADmax not specified. "<<endl;
        Info<<" Calculating from LADmax from LAI ... "<<endl;

        uniformLAD_ = false;
        LADmaxFromLAI();
    }

    calculateCanopy();
}

void Foam::fv::vegetationSource::correct()
{
    // solve the leaf heat balance at the end of each iteration
    // called in the top level solvers, if "hookOp solve;" is set.
    // Note that this is solved at the end of the iteration,
    // therefore, for iteration 0 Qs and Ql are 0
    leafHeatBalance();
}

Foam::tmp<Foam::volScalarField>
Foam::fv::vegetationSource::getRho() const
{

    if (obr_.foundObject<volScalarField>("rho"))
    {
        tmp<volScalarField> rho(obr_.lookupObject<volScalarField>("rho"));

        return rho;
    }
    else if (obr_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            obr_.lookupObject<transportModel>("transportProperties");

       tmp<volScalarField> rho(transport.rho());

       return rho;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get rho()" << nl
            << "Supported objects are: basicThermo  and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}

Foam::tmp<Foam::volScalarField>
Foam::fv::vegetationSource::getCp() const
{

    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo
            (obr_.lookupObject<fluidThermo>(basicThermo::dictName));

        tmp<volScalarField> Cp(thermo.Cp());

        return Cp;
    }
    else if (obr_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& transport =
            obr_.lookupObject<transportModel>("transportProperties");

        tmp<volScalarField> Cp(transport.Cp());

        return Cp;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find supported physical properties object to get Cp()" << nl
            << "Supported objects are: basicThermo and transportProperties"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


Foam::scalar Foam::fv::vegetationSource::setLAD(scalar z)
{
    if (z > height_ || z < 0.0)
    return 0;

    // calculate vertical distance represented by each value in LADProfile
    scalar deltaz = height_ / scalar(LADProfile_.size());

    // Calculate which index in the LADProfile that should be used at height z
    label profileIndex = mag(trunc(z / deltaz));

    // find index in LADProfile with max LAD
    label LADMaxIndex = 0;
    scalar tmpLAD = -1;
    for (label i=0; i < LADProfile_.size(); i++)
    {
        if (LADProfile_[i] > tmpLAD)
        {
            LADMaxIndex = i;
            tmpLAD = LADProfile_[LADMaxIndex];
        }
    }

    // return LAD at height z
    return LADProfile_[profileIndex] / LADProfile_[LADMaxIndex] * LADmax_;
}

void Foam::fv::vegetationSource::LADmaxFromLAI()
{
    LADmax_ = 0;
    scalar LAItmp = 0;
    while (LAItmp < LAI_)
    {
        LAItmp = integrateLAD();
        LADmax_ += 0.01;
    }
}

Foam::scalar Foam::fv::vegetationSource::integrateLAD()
{
    //uses simpsons rule to integrate between the points
    scalar LAItmp=0;
    for (label fi=1;fi<LADProfile_.size();fi++)
    {
        LAItmp += 0.5*(LADProfile_[fi-1]+LADProfile_[fi])*LADmax_*0.1*height_;
    }

    return LAItmp;
}

void Foam::fv::vegetationSource::leafHeatBalance()
{

    scalarField Tleaf(cells_.size(), 0);
    scalarField wleaf(cells_.size(), 0);

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);
    const volScalarField& p = obr_.lookupObject<volScalarField>(PName_);

    scalarField ra(cells_.size(), 0);   // aerodynamic resistance
    scalarField rs(cells_.size(), 150); // stomatal resistance

    // calculate aerodynamic resistance [s/m]
    forAll(cells_, i)
    {
        label celli = cells_[i];
        ra[i] = 130*pow((lsize_/max(mag(U[celli]), SMALL)), 0.5);
    }

    // get rhoCp fields
    tmp<volScalarField> rho(getRho());
    tmp<volScalarField> Cp(getCp());

    const volScalarField& G_  = obr_.lookupObject<volScalarField>("G");
    const volScalarField& Isolar_  = obr_.lookupObject<volScalarField>("Isolar");

    // Net absorbed radiative heat flux density in [w/m^3] = [1/m]*[W/m^2]
    scalarField Rn(cells_.size(), 0);

    forAll(cells_, i)
    {
        label celli = cells_[i];
        Rn[i]=(at_*G_[celli]+as_*Isolar_[celli]);
    }

    const volScalarField& w_  = obr_.lookupObject<volScalarField>("w");

    // absolute pressure in Pa
    scalarField pAb(cells_.size(), 0);

    if (obr_.foundObject<transportModel>("transportProperties"))
    {
        const IOdictionary& transportProperties
            = obr_.lookupObject<IOdictionary>("transportProperties");

        dimensionedScalar pRef("pRef",transportProperties.lookup("pRef"));
        scalar pRef_ = pRef.value();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            pAb[i] = p[celli]*rho()[celli] + pRef_; // Pa
        }
    }
    else if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const scalar pRef(obr_.lookupObject<fluidThermo>(basicThermo::dictName).pRefValue());
        forAll(cells_, i)
        {
            label celli = cells_[i];
            pAb[i] = p[celli] + pRef; // Pa
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot calculate the absolute pressure"
            << exit(FatalError);
    }

    // solve energy balance for each cell
    forAll(cells_, i)
    {
        label celli = cells_[i];

        // Newton's method
        scalar Tl_old = T[celli];
        scalar Tl_new = T[celli];

        scalar Tair = T[celli];

        label iloop = 0;
        scalar Tl_err = GREAT;
        scalar maxerr = 1e-6;
        label maxloops = 100;

        do
        {
            Tl_old = Tl_new;

            // hardcoded Tetens'equation for saturation pressure in Pa
            // requires temperature in degree C
            scalar Tl_old_c = Tl_old-273.15;

            scalar Psat = 610.78*Foam::exp((17.27*Tl_old_c)/(Tl_old_c+237.3));
            scalar dPsat = Psat*(17.27*(Tl_old_c+237.3)-17.27*Tl_old_c)/pow((Tl_old_c+237.3), 2);

            scalar wl = (Mvap_/Mair_)*Psat/pAb[i];

            scalar c1 = lambdaVap_*rho()[celli]/(ra[i]+rs[i]);
            scalar c2 = ra[i]/(2*rho()[celli]*Cp()[celli]);

            scalar Ql = c1*(wl-w_[celli]); // latent heat
            scalar dQl = c1*(Mvap_/Mair_)*dPsat/pAb[i];

            scalar f = Tair-Tl_old+c2*(Rn[i]/max(LAD[i], SMALL)-Ql);
            scalar df = -1 -c2*dQl;

            Tl_new = Tl_old - f/df;

            Tl_err = mag(Tl_new-Tl_old)/Tl_old;

            iloop++;
        }
        while (iloop < maxloops && Tl_err > maxerr);

        // update Tleaf
        Tleaf[i] = Tl_new;

        // update wleaf = f(Tleaf))
        scalar Tleaf_C = Tleaf[i]- 273.15;
        scalar Psat_new = 610.78*Foam::exp((17.27*Tleaf_C)/(Tleaf_C+237.3));
        wleaf[i] = (Mvap_/Mair_)*Psat_new/pAb[i];

        // Info<<" Ta "<<T[celli]<<" Tl "<<Tleaf[i]<<endl;
    }

    // calculate sensible and latent heat fluxes
    forAll(cells_, i)
    {
        label celli = cells_[i];
        Qs[i] = 2*rho()[celli]*Cp()[celli]*(Tleaf[i]-T[celli])/ra[i];
        Ql[i] = lambdaVap_*rho()[celli]*(wleaf[i]-w_[celli])/(ra[i]+rs[i]);
    }


}

void Foam::fv::vegetationSource::calculateCanopy()
{
    if (uniformLAD_)
    {
        forAll(cells_, i)
        {
            LAD[i] = LADmax_;
        }
    }
    else
    {
        scalarField distance(cells_.size(), 0);

        calculatePointDistance(distance);

        // set leaf area density profile for each cell
        // and store in class member LAD
        forAll(cells_, i)
        {
            scalar pointDistance = distance[i];
            if (pointDistance < height_)
            {
                LAD[i] = setLAD(pointDistance);
            }
        }
    }
    canopy_.reset
    (
       new volScalarField(
           IOobject
           (
              "canopy",
              mesh_.time().timeName(),
              obr_,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false  //not store in the register
           ),
           mesh_,
           dimensionedScalar("canopy",dimensionSet(0,-1,0,0,0,0,0),0)
       )
    );

    scalarField& canopy = canopy_->primitiveFieldRef();

    // calculates the canopy field as Cd_*LAD
    forAll(cells_, i)
    {
        label celli = cells_[i];
        canopy[celli] = Cd_*LAD[i];
    }
}

void Foam::fv::vegetationSource::calculatePointDistance(scalarField& distance)
{
    forAll(cells_, i)
    {
        // get global index
        label celli = cells_[i];
        const point cc = mesh_.cellCentres()[celli];
        vector v(cc - origin_);
        distance[i] = v & zoneDir_;
    }

    return;

}

// ************************************************************************* //

void Foam::fv::vegetationSource::addSup
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


void Foam::fv::vegetationSource::addSup
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

void Foam::fv::vegetationSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    // ico T
    if (eqn.psi().name() == word(TName_))
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

    // ico humidity
    if (eqn.psi().name() == word("w"))
    {
        const scalarField& V = mesh_.V();
        const volScalarField rho( getRho() );
        const volScalarField Cp( getCp() );

        scalarField& humiditySource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            humiditySource[celli] -=LAD[i]*Ql[i]/(rho[celli]*lambdaVap_)*V[celli];
        }
    }

}


void Foam::fv::vegetationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    // cmp energy eqn
    if (eqn.psi().name() == word(eName_))
    {
        const scalarField& V = mesh_.V();

        scalarField& heSource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            heSource[celli] -=LAD[i]*Qs[i]*V[celli];
        }
    }

    // cmp humidity
    if (eqn.psi().name() == word("w"))
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
