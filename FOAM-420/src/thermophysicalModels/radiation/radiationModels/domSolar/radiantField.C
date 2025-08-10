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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "radiationModels/domSolar/radiantField.H"
#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "cfdTools/general/bound/bound.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "cfdTools/general/include/fvCFD.H"
#include "cfdTools/general/explicitScalarWaveSolve/explicitScalarWaveSolve.H"
#include "fields/fvPatchFields/derived/radiant/radiantFvPatchScalarField.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "global/unitConversion/unitConversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* NamedEnum<radiantField::radiantInputType, 2>::names[] =
{
    "userDefined",
    "solarCalculator",
};

const NamedEnum<radiantField::radiantInputType, 2>
    radiantField::radiantInputTypeNames_(0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void radiantField::setFvDicts()
{
    word Iname = this->regIOobject::name();
    string pofix(string(Iname)(2, 100));

    fvMesh& refMesh = const_cast<fvMesh&>(this->mesh());

    // add solver options to relevant dictionaries as required
    dictionary fvSchemes = refMesh.schemes().localSchemeDict();

    // add div scheme
    {
        if (!fvSchemes.found(word("divSchemes")))
        {
            fvSchemes.add(word("divSchemes"), dictionary(), false);
        }

        char solarDiv[]
            = "bounded Gauss outletStabilised linearUpwind gradIenv";

        char divKeyc[] = "div(Jenv,Ienv)";

        fvSchemes.subDict("divSchemes").add
        (
            divKeyc,
            solarDiv,
            false
        );

        //add grad scheme for div scheme
        if (!fvSchemes.found(word("gradSchemes")))
        {
            fvSchemes.add(word("gradSchemes"), dictionary(), false);
        }

        char solarGrad[]
            = "cellLimited Gauss linear 1";

        char gradKeyc[] = "gradIenv";

        fvSchemes.subDict("gradSchemes").add
        (
            gradKeyc,
            solarGrad,
            false
        );


        refMesh.schemes().setLocalSchemeDict(fvSchemes);
    }

    dictionary fvSolution = refMesh.solution().localSolutionDict();

    // add solver settings
    {
        if (!fvSolution.found(word("solvers")))
        {
            fvSolution.add(word("solvers"), dictionary(), false);
        }

        dictionary Isolver;

        Isolver.set(word("solver"), word("PBiCG"));
        Isolver.set(word("preconditioner"), word("DILU"));
        Isolver.set(word("tolerance"), scalar(1e-05));
        Isolver.set(word("relTol"), scalar(0.1));

        fvSolution.subDict("solvers").add(word(Iname),Isolver);
        fvSolution.subDict("solvers").add(word(Iname+"Final"),Isolver);

        if (!fvSolution.found(word("relaxationFactors")))
        {
            fvSolution.add(word("relaxationFactors"), dictionary(), false);
            fvSolution.subDict("relaxationFactors").add(word("equations"), dictionary(), false);
        }
        else if (!fvSolution.subDict("relaxationFactors").found(word("equations")))
        {
            fvSolution.subDict("relaxationFactors").add(word("equations"), dictionary(), false);
        }

        dictionary& relaxFactor = fvSolution.subDict("relaxationFactors");
        relaxFactor.subDict("equations").add
        (
            Iname,
            scalar(1),
            false
        );
        refMesh.solution().setLocalSolutionDict(fvSolution);
    }

}

const dictionary radiantField::radiantBoundaryConditions
(
    const fvMesh& mesh,
    word fieldName
) const
{
    dictionary radiantBC;

    radiantBC.set(word("type"), word("radiant"));
    string pofix(string(fieldName)(1, 100));

    word JplusPofix = "J"+pofix;

    radiantBC.set(word("phi"), word("J"+pofix));
    char uvz[] = "uniform 0.0";
    radiantBC.set(word("inletValue"), uvz);
    radiantBC.set(word("value"), uvz);

    dictionary defaults;
    defaults.set(word("wall"), radiantBC);
    defaults.set(word("mappedWall"), radiantBC);
    defaults.set(word("indirectWall"), radiantBC);
    defaults.set(word("inlet"), radiantBC);
    defaults.set(word("outlet"), radiantBC);
    defaults.set(word("patch"), radiantBC);
    defaults.set(word("cyclicAMI"), radiantBC);

    dictionary bc;
    bc.add(word("boundaryTypeDefaults"), defaults, true);

    dictionary& subCyc
    (
        bc.subDict("boundaryTypeDefaults").subDict("cyclicAMI")
    );

    subCyc.add(word("patchType"), word("cyclicAMI"));

    return bc;
}

tmp<volScalarField> radiantField::autoCreateRadiantField
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    tmp<volScalarField> Ifield
    (
        new volScalarField
        (
            IOobject
            (
                "tmp"+fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(fieldName, dimMass/pow3(dimTime), 0.0)
        )
    );

    Ifield->boundaryFieldRef().reset
    (
        volScalarField::Boundary
        (
            mesh.boundary(),
            Ifield->internalField(),
            boundarySetup<scalar>
            (
                mesh,
                Ifield->internalField(),
                radiantBoundaryConditions(mesh, fieldName)
            )
        )
    );

    return Ifield;
}

scalar radiantField::inletOutletScaling()
{
    //calculate boundary flux so outlet flux can be scaled to match inlet flux
    const fvMesh& mesh = this->mesh();
    scalar inletFlux = 0;
    scalar outletFlux = 0;

    vector radiationDirection(0, 0, 0);

    switch (inputMode_)
    {
        case ritUserDefined :
        {
            radiationDirection = radiationDirection_->value(mesh.time().value());
            radiationDirection /= mag(radiationDirection);
        } break;

        case ritSolarCalculator :
        {
            radiationDirection = solarCalc_->direction();
        } break;
    }

    forAll(this->boundaryField(), bI)
    {
        const fvPatch& p = mesh.boundary()[bI];

        if (!p.coupled())
        {
            inletFlux += gSum
            (
                mag(neg(radiationDirection & p.nf())
                * ( radiationDirection & p.Sf() )
                * this->boundaryField()[bI])
            );
            outletFlux += gSum
            (
                mag(pos0(radiationDirection & p.nf())
                * ( radiationDirection & p.Sf() )
                * this->boundaryField()[bI])
            );
        }
    }

    //flux preserving scale
    scalar qsScale = Foam::min(Foam::max(inletFlux
        /Foam::max(SMALL, outletFlux), 0.9), 1.1);

    return qsScale;
}

void radiantField::updateRadiationProperties()
{
    //access mesh
    const fvMesh& mesh = this->mesh();


    switch (inputMode_)
    {
        case ritUserDefined :
        {
            currentFlux_ = radiationFlux_->value(mesh.time().value());
            currentDirection_
                = radiationDirection_->value(mesh.time().value());

        } break;

        case ritSolarCalculator :
        {
            switch (solarCalc_->sunDirectionModel())
            {
                case solarCalculator::mSunDirConstant:
                {
                    break;
                }
                case solarCalculator::mSunDirTracking:
                {
                    label updateIndex = label
                    (
                        mesh.time().value()
                        /solarCalc_->sunTrackingUpdateInterval()
                    );

                    if (updateIndex > updateTimeIndex_)
                    {
                        Info<< "Updating Sun position..." << endl;
                        updateTimeIndex_ = updateIndex;
                        solarCalc_->correctSunDirection();
                    }
                    break;
                }
            }

            currentFlux_ = solarCalc_->directSolarRad();
            currentDirection_ = solarCalc_->direction();

        Info<<"- Direct solar irradiation [W/m^2]: "<<currentFlux_<<endl;
        Info<<"- Solar Direction: "<<currentDirection_<<endl;

        } break;

        default:
        {
            FatalErrorInFunction
                << "Unsupported radiantInputType "
                << inputMode_
                << exit(FatalError);
        }
    }

    // set radiant field boundary conditions
    forAll(this->boundaryField(), bI)
    {
        if
        (
            isA<radiantFvPatchScalarField>
            (this->boundaryField()[bI])
        )
        {
            radiantFvPatchScalarField& cradp
            (
                dynamic_cast<radiantFvPatchScalarField&>
                (this->boundaryFieldRef()[bI])
            );

            cradp.radiationFlux() = currentFlux_;
            cradp.transmissivity() = patchSourceCoeffs_[bI];
        }
    }

    if (mag(currentDirection_) < VSMALL)
    {
        FatalErrorInFunction
            << "A zero sized vector has been specified for the direction of "
            << nl << "the radiant field " << internalField().name() << nl
            << exit(FatalError);
    }

    currentDirection_ /= mag(currentDirection_);
}

void radiantField::updateSkyDiffusiveRadiation
(
    volScalarField& qenv
)
{
    const fvMesh& mesh = this->mesh();

    vector verticalDir;

    switch (inputMode_)
    {
        case ritUserDefined :
        {
          // do nothing
        } break;

        case ritSolarCalculator :
        {
            switch (solarCalc_->sunLoadModel())
            {
                case solarCalculator::mSunLoadFairWeatherConditions:
                case solarCalculator::mSunLoadTheoreticalMaximum:
                {
                    verticalDir = verticalDirection();

                    const polyBoundaryMesh& patches = mesh.boundaryMesh();

                    forAll(patches, patchID)
                    {
                        const polyPatch& pp = patches[patchID];
                        if (isA<wallPolyPatch>(pp) || isA<mappedPatchBase>(pp))
                        {
                            const vectorField n = pp.faceNormals();

                            forAll(n, faceI)
                            {
                                const scalar cosEpsilon(verticalDir & -n[faceI]);

                                scalar Ed(0.0); // diffuse sky radiation
                                scalar Er(0.0); // diffuse ground reflected
                                const scalar cosTheta(solarCalc_->direction() & -n[faceI]);

                                // Above the horizon
                                if (cosEpsilon == 0.0)
                                {
                                    // Vertical walls
                                    scalar Y(0);

                                    if (cosTheta > -0.2)
                                    {
                                        Y = 0.55+0.437*cosTheta + 0.313*sqr(cosTheta);
                                    }
                                    else
                                    {
                                        Y = 0.45;
                                    }

                                    Ed = solarCalc_->C()*Y*solarCalc_->directSolarRad();
                                }
                                else
                                {
                                    // Other than vertical walls
                                    Ed =
                                        solarCalc_->C()
                                        * solarCalc_->directSolarRad()
                                        * (1.0 + cosEpsilon)/2.0;
                                }

                                // Ground reflected
                                Er =
                                    solarCalc_->directSolarRad()
                                    * (solarCalc_->C() + Foam::sin(solarCalc_->beta()))
                                    * solarCalc_->groundReflectivity()
                                    * (1.0 - cosEpsilon)/2.0;

                                qenv.boundaryFieldRef()[patchID][faceI] +=(Ed + Er);
                            }

                        } // if
                    }
                    break;
                }

                case solarCalculator::mSunLoadConstant:
                {
                    const polyBoundaryMesh& patches = mesh.boundaryMesh();

                    forAll(patches, patchID)
                    {
                        const polyPatch& pp = patches[patchID];
                        const vectorField n = pp.faceNormals();

                        if (isA<wallPolyPatch>(pp) || isA<mappedPatchBase>(pp))
                        {
                            forAll(n, faceI)
                            {
                                qenv.boundaryFieldRef()[patchID][faceI]
                                    +=solarCalc_->diffuseSolarRad();
                            }
                        }
                    }
                    break;
                }
            }
        } break;
    }
}

vector radiantField::verticalDirection()
{
    const fvMesh& mesh = this->mesh();

    vector vDir(0, 0, 0);

    switch (inputMode_)
    {
        case ritSolarCalculator :
        {
            vDir = solarCalc_->gridUp();
            vDir /= mag(vDir);
        } break;

        default:
        {
            if (mesh.foundObject<uniformDimensionedVectorField>("g"))
            {
                const uniformDimensionedVectorField& g =
                    mesh.lookupObject<uniformDimensionedVectorField>("g");
                vDir = (-g/mag(g)).value();
            }
            else
            {
                FatalErrorInFunction
                    << "Vertical direction gridUp not defined "
                    << exit(FatalError);
            }
        }
    }

    return vDir;
}


void radiantField::initialise(const dictionary& dict)
{
    //access mesh
    const fvMesh& mesh = this->mesh();

    //radiation input mode and components
    inputMode_ = radiantInputTypeNames_
        [dict.lookupOrDefault<word>("radiantInputType", "userDefined")];

    switch (inputMode_)
    {
        case ritUserDefined :
        {
            radiationFlux_.reset
            (
                Function1<scalar>::New("intensity", dict).ptr()
            );
            radiationDirection_.reset
            (
                Function1<vector>::New("direction", dict).ptr()
            );

        } break;

        case ritSolarCalculator :
        {
            solarCalc_.reset
            (
                new solarCalculator(dict, mesh)
            );

        } break;

        default:
        {
            FatalErrorInFunction
                << "Unsupported radiantInputType "
                << dict.lookupOrDefault<word>("radiantInputType", "userDefined")
                << exit(FatalError);
        }
    }


    updateTimeIndex_ = 0;
    initTol_ = dict.lookupOrDefault<scalar>("initTol", 1e-4);

    nSolverSolutions_ = dict.lookupOrDefault<label>("nSolutions", 20);
    solveField_ = dict.lookupOrDefault<bool>("solve", false);
    qScale_ = 1;

    IOobject IfieldHeader
    (
        this->name(),
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false //do not register
    );

    // read if present
    if (IfieldHeader.typeHeaderOk<volScalarField>(true))
    {
        this->forceAssign
        (
            tmp<volScalarField>
            (
                new volScalarField(IfieldHeader, mesh)
            )
        );

        //switches to determine initial behaviour based on state
        initializeField_ = false;
        updateAfterRead_ = true;
    }

    updateRadiationProperties();

    if (solveField_)
    {
        setFvDicts();
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

radiantField::radiantField
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const FieldField<Field, scalar>& sourcePatchCoeffs
)
:
    volScalarField
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            true
        ),
        autoCreateRadiantField
        (
            name,
            mesh
        )()
    ),

    inputMode_
    (
        radiantInputTypeNames_
        [dict.lookupOrDefault<word>("radiantInputType", "userDefined")]
    ),
    radiationFlux_(),
    radiationDirection_(),
    solarCalc_(),
    patchSourceCoeffs_(sourcePatchCoeffs),

    currentFlux_(0.0),
    currentDirection_(vector::zero),

    updateTimeIndex_(0),
    updateAfterRead_(false),
    initializeField_(true),

    initTol_
    (
        dict.lookupOrDefault<scalar>("initTol", 1e-4)
    ),
    nSolverSolutions_(dict.lookupOrDefault<label>("nSolutions", 20)),
    solveField_(dict.lookupOrDefault<bool>("solve", true)),
    qScale_(1)
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

radiantField::~radiantField(){}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void radiantField::calculate
(
    volScalarField& qenv
)
{
    const fvMesh& mesh = this->mesh();
    word Iname = this->regIOobject::name();
    string pofix(string(Iname)(1, 100));

    updateRadiationProperties();

    if (initializeField_) //only do this the first time or after update
    {
        surfaceScalarField Jflux
        (
            word("J"+pofix),
            (mesh.Sf() & currentDirection_)
        );

        SortableList<scalar> dist
        (
            (currentDirection_ & (mesh.C().internalField()))()
        );
        labelList order = dist.indices();

        scalar residual = 1;
        label iter = 0;

        do
        {
            this->correctBoundaryConditions();

            surfaceScalarField If( fvc::interpolate(*this) );

            Info<< iter++ << " ";
            residual = explicitScalarWaveSolve(*this, If, Jflux, order);

            reduce(residual, maxOp<scalar>());

        } while (residual > initTol_);

        qScale_ = inletOutletScaling();

        initializeField_ = false;
    }

    if (updateAfterRead_)
    {
        surfaceScalarField Jflux
        (
            word("J"+pofix),
            (mesh.Sf() & currentDirection_)
        );
        this->correctBoundaryConditions();
        qScale_ = inletOutletScaling();
        updateAfterRead_ = false;
    }

    if (solveField_ && mesh.time().timeIndex() < nSolverSolutions_)
    {
        // set radiant field boundary conditions
        forAll(this->boundaryField(), bI)
        {
            if
            (
                isA<radiantFvPatchScalarField>
                (this->boundaryField()[bI])
            )
            {
                radiantFvPatchScalarField& cradp
                (
                    dynamic_cast<radiantFvPatchScalarField&>
                    (this->boundaryFieldRef()[bI])
                );

                cradp.radiationFlux() = currentFlux_;
            }
        }


        surfaceScalarField Jflux
        (
            word("J"+pofix),
            (currentDirection_ & mesh.Sf() )
        );

        fvScalarMatrix IfieldEq
        (
            fvm::div(Jflux, *this, "div(Jenv,Ienv)")
        );

        //look up to see if there is a fvOptions in the registry
        if (mesh.thisDb().foundObject<fv::options>(fv::options::typeName))
        {
            fv::options& fvOptions = const_cast<fv::options&>
            (
            mesh.lookupObject<fv::options>(fv::options::typeName)
            );

            //dummy to pass dimensions check.
            tmp<volScalarField> dimfield
            (
            new volScalarField
            (
                IOobject
                (
                    "dimfield",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                mesh,
                dimensionedScalar("dimfield", dimTime/dimLength, 1.0)
             )
            );

            IfieldEq -= fvOptions(dimfield, *this);
            // IfieldEq -= fvOptions(*this);
        }

        IfieldEq.solve();

        bound
        (
            *this,
            dimensionedScalar(Iname, this->dimensions(), 0.0)
        );

        qScale_ = inletOutletScaling();

    }

    //- update boundary environmental radiant heat flux
    forAll(qenv.boundaryField(), bI)
    {
        const fvPatch& p = mesh.boundary()[bI];

        if (!p.coupled())
        {
            qenv.boundaryFieldRef()[bI]
                += pos0(currentDirection_ & p.nf())
                * (currentDirection_ & p.nf())
                * this->boundaryField()[bI]
                * qScale_;
        }
    }

    updateSkyDiffusiveRadiation(qenv);

}

Foam::scalar
Foam::radiantField::getZenithAngle() const
{
    switch (inputMode_)
        {
            case ritUserDefined :
            {
                const fvMesh& mesh = this->mesh();

                vector vDir(0, 0, 0);
                vector solarVectorDir(0, 0, 0);

                if (mesh.foundObject<uniformDimensionedVectorField>("g"))
                {
                    const uniformDimensionedVectorField& g =
                        mesh.lookupObject<uniformDimensionedVectorField>("g");
                    vDir = (-g/mag(g)).value();
                }
                else
                {
                    FatalErrorInFunction
                        << "gravity field g not found in order "
                        << " to determine vertical direction "
                        << exit(FatalError);
                }

                solarVectorDir = radiationDirection_->value(mesh.time().value());
                solarVectorDir /= mag(solarVectorDir);

                scalar cosPhi = (vDir & solarVectorDir)
                            /(mag(vDir)*mag(solarVectorDir)
                            + SMALL);

                // Enforce bounding between -1 and 1
                scalar radAngleBetween = acos(Foam::min(Foam::max(cosPhi, -1), 1));

                if (radAngleBetween > 1.57079632689) //90 degrees in rad
                {
                    radAngleBetween = 3.141592653 - radAngleBetween;
                }

                return radAngleBetween;

            } break;

            case ritSolarCalculator :
            {
                return 1.5707963 - solarCalc_->beta();
            } break;

            default:
            {
                FatalErrorInFunction
                    << "inputMode is invalid"
                    << exit(FatalError);
                return 0;
            }
        }

    return 0; //Dummy
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
