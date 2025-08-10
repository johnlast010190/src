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
    (c) 2010-2022 Esi Ltd.
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/porosityModel/DarcyForchheimer/DarcyForchheimer.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fvMatrices/fvMatrices.H"
#include "csv/CSVCore.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function1/xyzPolynomial/xyzPolynomial.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(DarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, DarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

Foam::scalar Foam::porosityModels::DarcyForchheimer::calcDcorrValue
(
    const scalar& T
) const
{
    return (aT1_ + aT2_*TRef_ + aT3_*sqr(TRef_))*(TRef_/T);
}


Foam::scalar Foam::porosityModels::DarcyForchheimer::calcFcorrValue
(
    const scalar& T
) const
{
    return (bT1_ + bT2_*TRef_ + bT3_*sqr(TRef_))*sqr(TRef_/T);
}


Foam::label Foam::porosityModels::DarcyForchheimer::fieldIndex(const label i) const
{
    label index = 0;
    if (!csys().uniform() || spatialDependancy_)
    {
        index = i;
    }
    return index;
}


void Foam::porosityModels::DarcyForchheimer::writeCoeffsFields() const
{
    volTensorField Dout
    (
        IOobject
        (
            typeName + ":D",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("0", dXYZ_.dimensions(), Zero)
    );
    volTensorField Fout
    (
        IOobject
        (
            typeName + ":F",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("0", fXYZ_.dimensions(), Zero)
    );

    forAll(cellZoneIDs_, zonei)
    {
        const label nCells = mesh_.cellZones()[cellZoneIDs_[zonei]].size();
        UIndirectList<tensor>(Dout, mesh_.cellZones()[cellZoneIDs_[zonei]]) =
            (D_[zonei].size() == 1) ? tensorField(nCells, D_[zonei].first()) : D_[zonei];
        UIndirectList<tensor>(Fout, mesh_.cellZones()[cellZoneIDs_[zonei]]) =
            (F_[zonei].size() == 1) ? tensorField(nCells, F_[zonei].first()) : F_[zonei];
    }

    Dout.write();
    Fout.write();
}


void Foam::porosityModels::DarcyForchheimer::initialise()
{
    if (coeffsMode_=="standard")
    {
        dimensionSet dDim = dimless/sqr(dimLength);
        dimensionSet fDim = dimless/dimLength;
        dXYZ_.dimensions().reset(dDim);
        fXYZ_.dimensions().reset(fDim);
        dXYZ_ = dimensionedVector("d", dDim, coeffs_);
        fXYZ_ = dimensionedVector("f", fDim, coeffs_);
    }
    else if (coeffsMode_=="pure")
    {
        //- currently alpha and beta are taken without dimensions.
        //  The reason for this is that they depend on the solver (p - p/rho)
        dXYZ_.name() = "alpha";
        fXYZ_.name() = "beta";

        if (pureCoeffsDef_=="userDefined")
        {
            // beta ~ kg/(m^3 s)
            //dimensionSet dDim = dimensionSet(1, -3, -1, 0, 0, 0, 0);
            // beta ~ kg/m^4
            //dimensionSet fDim = dimensionSet(1, -4, 0, 0, 0, 0, 0);
            dXYZ_.value() = vector(coeffs_.lookup(dXYZ_.name()));
            fXYZ_.value() = vector(coeffs_.lookup(fXYZ_.name()));
        }
        else if (pureCoeffsDef_=="cylindricalPolynomialFit")
        {
            if (csys().type()!="cylindrical")
            {
                //- check and throw error
                FatalErrorInFunction
                    << "cylindrical pure coeffs not compatible with "
                    << coorFramePtr_->coorSys().type()
                    << " coordinate system. "
                    << " Change it to cylidnrical."
                    << exit(FatalError);
            }
            if (pureCoeffsRhoMultiplier_)
            {
                FatalErrorInFunction
                    << "cylindrical pure coeffs not compatible with  "
                    << "pureCoeffsRhoMultiplier true."
                    << exit(FatalError);
            }
            if (spatialDependancy_)
            {
                FatalErrorInFunction
                    << "cylindricalPolynomialFit isn't compatible with "
                    << "spatialDependance true."
                    << exit(FatalError);
            }
            computeCylindricalCoeffs();
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported coefficients: "
            << "standard or pure are available"
            << exit(FatalError);
    }

    forAll(cellZoneIDs_, zonei)
    {
        const label zoneSize =
            mesh_.cellZones()[cellZoneIDs_[zonei]].size();

        Dcorr_[zonei].setSize(zoneSize, 1.0);
        Fcorr_[zonei].setSize(zoneSize, 1.0);
    }
    if (temperatureDependence_)
    {
        TName_ = coeffs_.lookupOrDefault<word>("T", "T");
        if (temperatureCoeffsMode_=="table")
        {
            ddn_ = Function1<scalar>::New("ddn", coeffs_);
            ffn_ = Function1<scalar>::New("ffn", coeffs_);
        }
        else if (temperatureCoeffsMode_=="polynomial")
        {
            setTemperatureDependenceCoeffs();
        }
        else
        {
            FatalErrorInFunction
                << " "
                << exit(FatalError);
        }
    }
    if (defaultMode_)
    {
        adjustNegativeResistance(dXYZ_);
        adjustNegativeResistance(fXYZ_);
    }
    if (spatialDependancy_)
    {
        dSpatialFun_ = Function1<vector>::New("dSpatialCoeffs", coeffs_);
        fSpatialFun_ = Function1<vector>::New("fSpatialCoeffs", coeffs_);
        const word xyzTypeName(Function1Types::xyzPolynomial<vector>::typeName);
        if
        (
            dSpatialFun_->type() != xyzTypeName
         || fSpatialFun_->type() != xyzTypeName
        )
        {
            FatalErrorInFunction
                << Foam::porosityModels::DarcyForchheimer::type()
                << " for the has unsupported spatial fuction name: \""
                << (
                      (dSpatialFun_->type() != xyzTypeName)
                    ? dSpatialFun_->type()
                    : fSpatialFun_->type()
                  ) << "\".\n"
                << "The only supported type is \"" << xyzTypeName << "\""
                << exit(FatalError);
        }

        updateSpatialCoeffs();
    }
    else
    {
        calcTransformModelData();
    }
}


void Foam::porosityModels::DarcyForchheimer::computeCylindricalCoeffs()
{
    fileName fName(coeffs_.lookup("file"));
    if (fName.ext()!= "csv")
    {
        FatalErrorInFunction
            << "Only csv format is supported for dp and Q file definition. "
            << exit(FatalError);
    }

    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fName.expand()));
    ISstream& is = isPtr();

    if (!is.good())
    {
        FatalIOErrorInFunction(is) << "Problem with CSV file."
            << exit(FatalIOError);
    }


    Switch mergeSeparators =
        coeffs_.lookupOrDefault<Switch>("mergeSeparators", false);
    char separator =
        coeffs_.lookupOrDefault<string>("separator", string(","))[0];

    //- Header read
    string headerLine;
    is.getLine(headerLine);

    List<string> header
    (
        fileFormats::CSVCore::readLine(headerLine, separator, mergeSeparators)
    );
    labelList loc(2, label(-1));
    forAll(header, sI)
    {
        if (header[sI] == "dp") loc[0] = sI;
        if (header[sI] == "Q") loc[1] = sI;
    }

    DynamicList<scalar> dpSample;
    DynamicList<scalar> qSample;
    while (is.good())
    {
        string line;
        is.getLine(line);

        List<string> splitted
        (
            fileFormats::CSVCore::readLine(line, separator, mergeSeparators)
        );
        if (splitted.size() <= 1) break;

        //- coordinate
        {
            scalar dpI = readScalar(IStringStream(splitted[loc[0]])());
            scalar qI = readScalar(IStringStream(splitted[loc[1]])());
            dpSample.append(dpI);
            qSample.append(qI);
        }
    }
    dpSample.shrink();
    qSample.shrink();

    dpList_.transfer(dpSample);
    qList_.transfer(qSample);

    scalar rad1 = readScalar(coeffs_.lookup("innerRadius"));
    scalar rad2 = readScalar(coeffs_.lookup("outerRadius"));
    scalar height = readScalar(coeffs_.lookup("height"));
    scalar sectors = coeffs_.lookupOrDefault<scalar>("sectors", 1);
    scalar rhoRef = readScalar(coeffs_.lookup("rhoRef"));

    scalar C1 = 0.0;
    scalar C2 = 0.0;

    Info<< endl;
    Info<< "Polynomial fit with dp:" << endl;
    polynomialFit(C1, C2, rhoRef, rad1, rad2, height, sectors);

    Info<< "dp(U) = " << -C1 << " U " << -C2 << " U^2"  << endl;
    Info<< endl;

    dXYZ_.value() = vector(C1, 0.0, 0.0);
    fXYZ_.value() = vector(C2, 0.0, 0.0);
}


void Foam::porosityModels::DarcyForchheimer::polynomialFit
(
    scalar& C1,
    scalar& C2,
    const scalar& rhoRef,
    const scalar& r1,
    const scalar& r2,
    const scalar& height,
    const scalar& sectors
) const
{

    const scalar pi = constant::mathematical::pi;

    //Polynomial regression y = a1*x+a2*x^2
    //creating the matrix to solve the system
    // sum[x] sum[x^2] sum[x^3] sum[x^4] sum[YX] sum[y x^2]
    scalar sumXp2 = 0.0;
    scalar sumXp3 = 0.0;
    scalar sumXp4 = 0.0;
    scalar sumYX = 0.0;
    scalar sumYXp2 = 0.0;
    forAll(qList_, pI)
    {
        scalar xp2I = sqr(qList_[pI]);
        sumXp2 += xp2I;
        sumXp3 += xp2I*qList_[pI];
        sumXp4 += sqr(xp2I);
        sumYX += dpList_[pI]*qList_[pI];
        sumYXp2 += dpList_[pI]*xp2I;
    }

    if (debug)
    {
        Info<< "coeffs: " << tab
             << sumXp2 << tab
             << sumXp3 << tab
             << sumXp4 << tab
             << sumYX << tab
             << sumYXp2 << tab
             << endl;
    }

    scalar a2 = (sumYX-sumYXp2*sumXp2/sumXp3)/(sumXp3-sumXp2/sumXp3);
    scalar a1 = (sumYX-sumXp3*a2)/(sumXp2);

    Info<< "dp(Q) = " << a1 << " Q " << a2 << " Q^2"  << endl;

    scalar aMul1 = 1/(rhoRef*2.0*pi*height/sectors);
    scalar aMul2 = sqr(aMul1);

    scalar aMul3 = log(r2/r1);
    scalar aMul4 = 1/r1 - 1/r2;

    C1 = -a1/(aMul1*aMul3);
    C2 = -a2/(aMul2*aMul4);

    if (debug)
    {
        Info<< "C1: " << C1 << tab << "C2: " << C2 << endl;
    }
}


void Foam::porosityModels::DarcyForchheimer::updateCorrCoeffs() const
{
    bool updated = false;
    if (spatialDependancy_ || temperatureDependence_)
    {
        updated = (curTimeIndex_== mesh_.time().timeIndex());
    }
    // Update spatial coefficients if mesh is moving
    if (spatialDependancy_ && !updated && mesh_.changing())
    {
        updateSpatialCoeffs();
    }
    if (temperatureDependence_ && obr_.foundObject<volScalarField>(TName_) && !updated)
    {
        const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

        forAll(cellZoneIDs_, zoneI)
        {
            scalarField& ddnZones = Dcorr_[zoneI];
            scalarField& ffnZones = Fcorr_[zoneI];
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
            if (temperatureCoeffsMode_=="table")
            {
                forAll(cells, i)
                {
                    label cellI = cells[i];
                    ddnZones[i] = ddn_->value(T[cellI]);
                    ffnZones[i] = ffn_->value(T[cellI]);
                }
            }
            else
            {
                forAll(cells, i)
                {
                    label cellI = cells[i];
                    ddnZones[i] = calcDcorrValue(T[cellI]);
                    ffnZones[i] = calcFcorrValue(T[cellI]);
                }
            }
        }
    }

    if (debug && mesh_.time().writeTime() && !updated)
    {
        writeCoeffsFields();
    }

    if
    (
        (
            spatialDependancy_
         || (temperatureDependence_ && obr_.foundObject<volScalarField>(TName_))
        )
     && !updated)
    {
        curTimeIndex_ = mesh_.time().timeIndex();
    }
}


void Foam::porosityModels::DarcyForchheimer::setTemperatureDependenceCoeffs()
{
    TRef_ = coeffs_.lookupOrDefault<scalar>("TRef", 0.0);
    aT1_ = coeffs_.lookupOrDefault<scalar>("a1", 0.0);
    aT2_ = coeffs_.lookupOrDefault<scalar>("a2", 0.0);
    aT3_ = coeffs_.lookupOrDefault<scalar>("a3", 0.0);
    bT1_ = coeffs_.lookupOrDefault<scalar>("b1", 0.0);
    bT2_ = coeffs_.lookupOrDefault<scalar>("b2", 0.0);
    bT3_ = coeffs_.lookupOrDefault<scalar>("b3", 0.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::DarcyForchheimer
(
    const word& name,
    const word& modelType,
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, obr, mesh, dict, cellZoneName),
    dXYZ_("tmpD", dimless, vector::zero),
    fXYZ_("tmpF", dimless, vector::zero),
    D_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    Dcorr_(cellZoneIDs_.size()),
    Fcorr_(cellZoneIDs_.size()),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "nu")),
    defaultMode_(coeffs_.lookupOrDefault<Switch>("defaultMode", true)),
    coeffsMode_(coeffs_.lookupOrDefault<word>("coeffsMode", "standard")),
    pureCoeffsDef_
    (
        coeffs_.lookupOrDefault<word>("pureCoeffsDefinition", "userDefined")
    ),
    pureCoeffsRhoMultiplier_
    (
        coeffs_.lookupOrDefault<Switch>("pureCoeffsRhoMultiplier", "false")
    ),
    temperatureDependence_
    (
        coeffs_.lookupOrDefault<Switch>("temperatureDependence", "false")
    ),
    temperatureCoeffsMode_
    (
        coeffs_.lookupOrDefault<word>("temperatureCoeffsMode", "polynomial")
    ),
    spatialDependancy_
    (
        coeffs_.lookupOrDefault<Switch>("spatialDependance", false)
    ),
    TName_(word::null),
    ddn_(),
    ffn_(),
    TRef_(0.0),
    aT1_(0.0),
    aT2_(0.0),
    aT3_(0.0),
    bT1_(0.0),
    bT2_(0.0),
    bT3_(0.0),
    dpList_(0),
    qList_(0),
    curTimeIndex_(-1)
{
    if (modelType == "temperatureDependentDarcyForchheimer")
    {
        FatalErrorInFunction
            << "temperatureDependentDarcyForchheimer "
            << "has been merged with DarcyForchheimer class. "
            << "Check DarcyForchheimer different modes."
            << exit(FatalError);
    }
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::~DarcyForchheimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::DarcyForchheimer::updateSpatialCoeffs() const
{
    const labelList diagonal({tensor::XX, tensor::YY, tensor::ZZ});
    tensor darcyCoeff(Zero);
    tensor forchCoeff(Zero);
    forAll(diagonal, i)
    {
        const label compI = diagonal[i];
        darcyCoeff[compI] = dXYZ_.value()[i];
        forchCoeff[compI] = fXYZ_.value()[i];
        if (coeffsMode_ == "standard")
        {
            forchCoeff[compI] *= 0.5;
        }
    }
    forAll(cellZoneIDs_, zonei)
    {
        const pointField cc
        (
            mesh_.cellCentres(),
            mesh_.cellZones()[cellZoneIDs_[zonei]]
        );

        const vectorField ccLocal(csys().localPosition(cc));
        const vectorField dCoeffs
        (
            csys().globalPosition(dSpatialFun_->value(ccLocal)())
        );
        const vectorField fCoeffs
        (
            csys().globalPosition(fSpatialFun_->value(ccLocal)())
        );

        // Replacement for xx, yy, zz works only in
        // global coordinate system
        tensorField tdCoeffs(dCoeffs.size(), Zero);
        tensorField tfCoeffs(fCoeffs.size(), Zero);
        forAll(diagonal, i)
        {
            const label compI = diagonal[i];
            tdCoeffs.replace(compI, dCoeffs.component(i));
            tfCoeffs.replace(compI, fCoeffs.component(i));
        }

        D_[zonei] = csys().transform(cc, cmptMultiply(darcyCoeff, tdCoeffs));
        F_[zonei] = csys().transform(cc, cmptMultiply(forchCoeff, tfCoeffs));
    }
}


void Foam::porosityModels::DarcyForchheimer::calcTransformModelData()
{
    // The Darcy coefficient as a tensor
    tensor darcyCoeff(Zero);
    darcyCoeff.xx() = dXYZ_.value().x();
    darcyCoeff.yy() = dXYZ_.value().y();
    darcyCoeff.zz() = dXYZ_.value().z();

    // The Forchheimer coefficient as a tensor
    tensor forchCoeff(Zero);
    forchCoeff.xx() = fXYZ_.value().x();
    forchCoeff.yy() = fXYZ_.value().y();
    forchCoeff.zz() = fXYZ_.value().z();

    if (coeffsMode_=="standard")
    {
        // - the leading 0.5 is from 1/2*rho
        forchCoeff.xx() *= 0.5;
        forchCoeff.yy() *= 0.5;
        forchCoeff.zz() *= 0.5;
    }

    if (csys().uniform())
    {
        forAll(cellZoneIDs_, zonei)
        {
            D_[zonei].resize(1);
            F_[zonei].resize(1);

            D_[zonei] = csys().transform(darcyCoeff);
            F_[zonei] = csys().transform(forchCoeff);
        }
    }
    else
    {
        forAll(cellZoneIDs_, zonei)
        {
            const pointUIndList cc
            (
                mesh_.cellCentres(),
                mesh_.cellZones()[cellZoneIDs_[zonei]]
            );
            D_[zonei] = csys().transform(cc, darcyCoeff);
            F_[zonei] = csys().transform(cc, forchCoeff);
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    //- add it for postProcess mode
    updateCorrCoeffs();

    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    if (coeffsMode_=="standard")
    {
        apply(Udiag, Usource, V, rho, mu, U);
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {
            apply(Udiag, Usource, V, rho, rho, U);
        }
        else
        {
            scalarField oneField(Udiag.size(), 1.0);
            apply(Udiag, Usource, V, geometricOneField(), oneField, U);
        }
    }

    force = Udiag*U - Usource;
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    updateCorrCoeffs();

    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    if (coeffsMode_=="standard")
    {
        word rhoName(IOobject::groupName(rhoName_, U.group()));
        word muName(IOobject::groupName(muName_, U.group()));
        word nuName(IOobject::groupName(nuName_, U.group()));

        if (UEqn.dimensions() == dimForce)
        {
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>(rhoName);

            if (obr_.foundObject<volScalarField>(muName))
            {
                const volScalarField& mu =
                    obr_.lookupObject<volScalarField>(muName);

                apply(Udiag, Usource, V, rho, mu, U);
            }
            else
            {
                const volScalarField& nu =
                    obr_.lookupObject<volScalarField>(nuName);

                apply(Udiag, Usource, V, rho, rho*nu, U);
            }
        }
        else
        {
            if (obr_.foundObject<volScalarField>(nuName))
            {
                const volScalarField& nu =
                    obr_.lookupObject<volScalarField>(nuName);

                apply(Udiag, Usource, V, geometricOneField(), nu, U);
            }
            else
            {
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);
                const volScalarField& mu =
                    obr_.lookupObject<volScalarField>(muName);

                apply(Udiag, Usource, V, geometricOneField(), mu/rho, U);
            }
        }
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {

            if (UEqn.dimensions() == dimForce)
            {
                word rhoName(IOobject::groupName(rhoName_, U.group()));

                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                apply(Udiag, Usource, V, rho, rho, U);
            }
            else
            {
                // ico source terms is dP/(dx*rho)
                scalarField oneField(Udiag.size(), 1.0);
                apply(Udiag, Usource, V, geometricOneField(), oneField, U);
            }
        }
        else
        {
            scalarField oneField(Udiag.size(), 1.0);

            if (UEqn.dimensions() == dimForce)
            {
                apply(Udiag, Usource, V, geometricOneField(), oneField, U);
            }
            else
            {
                word rhoName(IOobject::groupName(rhoName_, U.group()));

                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                // ico source terms is dP/(dx*rho)
                apply(Udiag, Usource, V, (1/rho)(), (1/rho)(), U);
            }
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvBlockMatrix<vector>& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    vectorField& Usource = UEqn.source();

    updateCorrCoeffs();

    if (coeffsMode_=="standard")
    {
        word rhoName(IOobject::groupName(rhoName_, U.group()));
        word muName(IOobject::groupName(muName_, U.group()));
        word nuName(IOobject::groupName(nuName_, U.group()));

        if (UEqn.dimensionSets()[0] == dimForce)
        {
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>(rhoName);
            if (obr_.foundObject<volScalarField>(nuName))
            {
                const volScalarField& nu =
                    obr_.lookupObject<volScalarField>(nuName);

                apply(blockDiag, Usource, rho, rho*nu, U);
            }
            else
            {
                const volScalarField& mu =
                    obr_.lookupObject<volScalarField>(muName);

                apply(blockDiag, Usource, rho, mu, U);
            }
        }
        else
        {
            if (obr_.foundObject<volScalarField>(nuName))
            {
                const volScalarField& nu =
                    obr_.lookupObject<volScalarField>(nuName);

                apply(blockDiag, Usource, geometricOneField(), nu, U);
            }
            else
            {
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);
                const volScalarField& mu =
                    obr_.lookupObject<volScalarField>(muName);

                apply(blockDiag, Usource, geometricOneField(), mu/rho, U);
            }
        }
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {
            scalarField oneField(Usource.size(), 1.0);

            if (UEqn.dimensionSets()[0] == dimForce)
            {
                word rhoName(IOobject::groupName(rhoName_, U.group()));
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                apply(blockDiag, Usource, rho, rho, U);
            }
            else
            {
                apply(blockDiag, Usource, geometricOneField(), oneField , U);
            }
        }
        else
        {
            scalarField oneField(Usource.size(), 1.0);

            if (UEqn.dimensionSets()[0] == dimForce)
            {
                apply(blockDiag, Usource, geometricOneField(), oneField, U);
            }
            else
            {
                word rhoName(IOobject::groupName(rhoName_, U.group()));

                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                // ico source terms is dP/(dx*rho)
                apply(blockDiag, Usource, (1/rho)(), (1/rho)(), U);
            }
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    updateCorrCoeffs();

    if (coeffsMode_=="standard")
    {
        apply(Udiag, Usource, V, rho, mu, U);
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {
            scalarField oneField(Udiag.size(), 1.0);
            apply(Udiag, Usource, V, rho, rho, U);
        }
        else
        {
            scalarField oneField(Udiag.size(), 1.0);
            apply(Udiag, Usource, V, geometricOneField(), oneField, U);
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvBlockMatrix<vector>& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const volVectorField& U = UEqn.psi();

    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    vectorField& Usource = UEqn.source();

    updateCorrCoeffs();

    if (coeffsMode_=="standard")
    {
        apply(blockDiag, Usource, rho, mu, U);
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {
            scalarField oneField(Usource.size(), 1.0);
            apply(blockDiag, Usource, rho, rho, U);
        }
        else
        {
            scalarField oneField(Usource.size(), 1.0);
            apply(blockDiag, Usource, geometricOneField(), oneField, U);
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    updateCorrCoeffs();
    vectorField nullSource (AU.size(), 0);

    if (coeffsMode_=="standard")
    {
        word rhoName(IOobject::groupName(rhoName_, U.group()));
        word muName(IOobject::groupName(muName_, U.group()));
        word nuName(IOobject::groupName(nuName_, U.group()));

        if (UEqn.dimensions() == dimForce)
        {
            const volScalarField& rho
                = obr_.lookupObject<volScalarField>(rhoName);

            if (obr_.foundObject<volScalarField>(muName))
            {
                const volScalarField& mu
                    = obr_.lookupObject<volScalarField>(muName);

                apply(AU, nullSource, rho, mu, U);
            }
            else
            {
                const volScalarField& nu
                    = obr_.lookupObject<volScalarField>(nuName);

                apply(AU, nullSource, rho, rho*nu, U);
            }
        }
        else
        {
            if (obr_.foundObject<volScalarField>(nuName))
            {
                const volScalarField& nu =
                    obr_.lookupObject<volScalarField>(nuName);

                apply(AU, nullSource, geometricOneField(), nu, U);
            }
            else
            {
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);
                const volScalarField& mu =
                    obr_.lookupObject<volScalarField>(muName);

                apply(AU, nullSource, geometricOneField(), mu/rho, U);
            }
        }
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {

            if (UEqn.dimensions() == dimForce)
            {
                word rhoName(IOobject::groupName(rhoName_, U.group()));
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                apply(AU, nullSource, rho, rho, U);
            }
            else
            {
                scalarField oneField(UEqn.diag().size(), 1.0);
                apply(AU, nullSource, geometricOneField(), oneField, U);
            }
        }
        else
        {
            scalarField oneField(UEqn.diag().size(), 1.0);

            if (UEqn.dimensions() == dimForce)
            {
                apply(AU, nullSource, geometricOneField(), oneField, U);
            }
            else
            {
                word rhoName(IOobject::groupName(rhoName_, U.group()));

                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                apply(AU, nullSource, (1/rho)(), (1/rho)(), U);
            }
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volVectorField& Uprimal
) const
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    updateCorrCoeffs();

    const scalarField& V = mesh_.V();
    scalarField& Udiag = UaEqn.diag();
    vectorField& Usource = UaEqn.source();

    dimensionSet idDim
        = UaEqn.dimensions()*dimTime/dimVolume/UaEqn.psi().dimensions();

    const volVectorField& U = UaEqn.psi();

    if (coeffsMode_=="standard")
    {
        word rhoName(IOobject::groupName(rhoName_, Uprimal.group()));
        word muName(IOobject::groupName(muName_, Uprimal.group()));
        word nuName(IOobject::groupName(nuName_, Uprimal.group()));

        if (idDim == dimDensity)
        {
            const volScalarField& rho
                = obr_.lookupObject<volScalarField>(rhoName_);

            if (obr_.foundObject<volScalarField>(muName))
            {
                const volScalarField& mu
                    = obr_.lookupObject<volScalarField>(muName);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    rho,
                    mu,
                    Uprimal,
                    U
                );
            }
            else
            {
                const volScalarField& nu
                    = obr_.lookupObject<volScalarField>(nuName);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    rho,
                    rho*nu,
                    Uprimal,
                    U
                );
            }
        }
        else if (idDim == dimless)
        {
            if (obr_.foundObject<volScalarField>(nuName))
            {
                const volScalarField& nu =
                    obr_.lookupObject<volScalarField>(nuName);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    geometricOneField(),
                    nu,
                    Uprimal,
                    U
                );
            }
            else
            {
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);
                const volScalarField& mu =
                    obr_.lookupObject<volScalarField>(muName);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    geometricOneField(),
                    mu/rho,
                    Uprimal,
                    U
                );
            }
        }
        else
        {
            //unsupported dimensions
            FatalErrorInFunction
                << "Unsupported adjoint matrix dimensions: "
                << UaEqn.dimensions() << exit(FatalError);
        }
    }
    else
    {
        if (pureCoeffsRhoMultiplier_)
        {
            if (idDim == dimDensity)
            {
                word rhoName(IOobject::groupName(rhoName_, Uprimal.group()));
                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    rho,
                    rho,
                    Uprimal,
                    U
                );
            }
            else
            {

                scalarField oneField(Usource.size(), 1.0);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    geometricOneField(),
                    oneField,
                    Uprimal,
                    U
                );
            }
        }
        else
        {
            if (idDim == dimDensity)
            {
                scalarField oneField(Usource.size(), 1.0);
                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    geometricOneField(),
                    oneField,
                    Uprimal,
                    U
                );
            }
            else
            {
                word rhoName(IOobject::groupName(rhoName_, Uprimal.group()));

                const volScalarField& rho =
                    obr_.lookupObject<volScalarField>(rhoName);

                adjointApply
                (
                    Udiag,
                    Usource,
                    V,
                    (1/rho)(),
                    (1/rho)(),
                    Uprimal,
                    U
                );
            }
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Uprimal
) const
{
    updateCorrCoeffs();

    adjointCorrect(UaEqn, Uprimal);
}


void Foam::porosityModels::DarcyForchheimer::adjointCorrect
(
    fvBlockMatrix<vector>& UaEqn,
    const volVectorField& U
) const
{
    updateCorrCoeffs();

    const scalarField& V = mesh_.V();

    typename CoeffField<vector>::squareTypeField& diag =
        UaEqn.diag().asSquare();

    if (coeffsMode_=="standard")
    {
        const volScalarField& nu = obr_.lookupObject<volScalarField>(nuName_);

        forAll(cellZoneIDs_, zoneI)
        {
            tensorField dZonest = D_[zoneI].T()();
            tensorField fZonest = F_[zoneI].T()();
            const tensorField& fZones = F_[zoneI];

            const scalarField& ddnZones = Dcorr_[zoneI];
            const scalarField& ffnZones = Fcorr_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                const label j = fieldIndex(i);

                vector Upr(U[cellI]);
                if (coorFramePtr_)
                {
                    Upr -= coorFramePtr_->frameVelocity(mesh_.C()[cellI], true);
                }
                scalar magUprimal = max(mag(Upr), SMALL);

                const tensor adjDragCoeff =
                    max(
                        (
                            nu[cellI]*dZonest[j]*ddnZones[i]
                          + (magUprimal)*fZonest[j]*ffnZones[i]
                          + (
                                1/(magUprimal)*(fZones[j]*
                                ffnZones[i] & sqr(Upr)).T()
                            )
                        )*V[cellI],
                        pTraits<tensor>::zero
                    );

                diag[cellI](0, 0) += adjDragCoeff.xx();
                diag[cellI](0, 1) += adjDragCoeff.xy();
                diag[cellI](0, 2) += adjDragCoeff.xz();
                diag[cellI](1, 0) += adjDragCoeff.yx();
                diag[cellI](1, 1) += adjDragCoeff.yy();
                diag[cellI](1, 2) += adjDragCoeff.yz();
                diag[cellI](2, 0) += adjDragCoeff.zx();
                diag[cellI](2, 1) += adjDragCoeff.zy();
                diag[cellI](2, 2) += adjDragCoeff.zz();
            }
        }
    }
    else
    {

        forAll(cellZoneIDs_, zoneI)
        {
            tensorField dZonest = D_[zoneI].T()();
            tensorField fZonest = F_[zoneI].T()();
            const tensorField& fZones = F_[zoneI];

            const scalarField& ddnZones = Dcorr_[zoneI];
            const scalarField& ffnZones = Fcorr_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            const volScalarField& rho =
                obr_.lookupObject<volScalarField>(rhoName_);

            forAll(cells, i)
            {
                const label cellI = cells[i];
                const label j = fieldIndex(i);

                vector Upr(U[cellI]);
                if (coorFramePtr_)
                {
                    Upr -= coorFramePtr_->frameVelocity(mesh_.C()[cellI], true);
                }
                scalar magUprimal = max(mag(Upr), SMALL);

                const tensor adjDragCoeff =
                    (
                          dZonest[j]*ddnZones[i]/rho[cellI]
                        + (magUprimal)*fZonest[j]*ffnZones[i]/rho[cellI]
                        + (
                            1/(magUprimal)/rho[cellI]*
                            (fZones[j]*ffnZones[i] & sqr(Upr)).T()
                          )
                    )*V[cellI];

                diag[cellI](0, 0) += adjDragCoeff.xx();
                diag[cellI](0, 1) += adjDragCoeff.xy();
                diag[cellI](0, 2) += adjDragCoeff.xz();
                diag[cellI](1, 0) += adjDragCoeff.yx();
                diag[cellI](1, 1) += adjDragCoeff.yy();
                diag[cellI](1, 2) += adjDragCoeff.yz();
                diag[cellI](2, 0) += adjDragCoeff.zx();
                diag[cellI](2, 1) += adjDragCoeff.zy();
                diag[cellI](2, 2) += adjDragCoeff.zz();
            }
        }
    }
}


bool Foam::porosityModels::DarcyForchheimer::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
