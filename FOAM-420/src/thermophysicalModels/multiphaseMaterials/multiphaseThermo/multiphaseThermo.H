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
    (c) 2023 Esi Ltd.

Class
    Foam::multiphaseThermo

SourceFiles
    multiphaseThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseThermo_H
#define multiphaseThermo_H

#include "fields/volFields/volFields.H"
#include "psiThermo/psiThermo.H"
#include "rhoThermo/rhoThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "mixtures/phaseVolumeFractions/phaseVolumeFractions.H"
#include "interfacePair/interfacePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class multiphaseThermo Declaration
\*---------------------------------------------------------------------------*/

class multiphaseThermo
:
    public rhoThermo
{
    // Private Data

        //- Reference to material models
        materialTables& materials_;

        //- Pointer list for phases thermos
        PtrList<basicThermo> thermos_;

        //- Phase list
        const wordList phases_;

        //- Pointer to volume fraction class
        phaseVolumeFractions* fractions_;


    // Private Member Functions

        //- Construct material tables
        materialTables& matLookupOrConstruct
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Return boundary field for different combination of
        // temperature/pressure
        tmp<scalarField> tmpBoundaryField
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi,
            const word& modelName
        ) const;


public:

    //- Runtime type information
    TypeName("multiphase");


    // Constructors

        //- Construct from object registry
        // (phase name not applicable but necessary
        // to add to run time selection tables)
        multiphaseThermo
        (
            const objectRegistry& obr,
            const word& phaseName = word::null
        );


    // Selectors

        //- Generic lookup or create for any thermo type
        //  Returns the sub-thermo for the specified phaseName if given
        static basicThermo& lookupOrCreate
        (
            const objectRegistry&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~multiphaseThermo()
    {}


    // Member Functions

        //- Additional initialisation - called by parent heThermo after
        //  construction
        virtual void init();

        //- Update properties
        virtual void correct();

        //- Update densities for given pressure change
        void correctRho(const volScalarField& dp);

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;

        //- Thermos for phases
        const PtrList<basicThermo>& thermos() const
        {
            return thermos_;
        }

        //- Reference to thermos for phases
        PtrList<basicThermo>& thermos()
        {
            return thermos_;
        }

        //- Thermos for phases
        const PtrList<volScalarField>& alphas() const
        {
            return fractions_->alphas();
        }

        //- Thermos for phases
        PtrList<volScalarField>& alphas()
        {
            return fractions_->alphas();
        }

        //- Thermos for pahses
        PtrList<volScalarField>& dgdts() const
        {
            return const_cast<PtrList<volScalarField>&>(fractions_->dgdts());
        }

        const wordList& phases() const
        {
            return phases_;
        }

        phaseVolumeFractions& fractions()
        {
            FOAM_ASSERT(fractions_)
            {
                FatalErrorInFunction
                    << "Volume fractions is uninitialised"
                    << exit(FatalError);
            }
            return *fractions_;
        }

        // Fields derived from thermodynamic state variables

            //- Return the phase-averaged reciprocal Cv
            tmp<volScalarField> rCv() const;


    // TODO Functions and data to be moved

        //- Get mesh through object/solution registry
        const fvMesh& getMesh(const objectRegistry& obr) const;

        //- Mesh reference.
        //TODO: should be probably taken from getMesh func?
        const fvMesh& mesh_;

        //- Hash table for sigma interphase pairs
        typedef HashPtrTable<Function1<scalar>, interfacePair, interfacePair::hash>
            sigmaTable;

        //- Helper class to create Function1's from stream for HashPtrTable
        class iNewSigma
        {
        public:

            iNewSigma()
            {}

            autoPtr<Function1<scalar>> operator()
            (
                const interfacePair& key, Istream& is
            ) const
            {
                word keyName = "sigma";
                token firstToken(is);
                if (firstToken.isScalar())
                {
                    dictionary d;
                    d.add(keyName, firstToken);
                    return
                        autoPtr<Function1<scalar>>
                        (
                            Function1<scalar>::New(keyName, d)
                        );
                }
                is.putBack(firstToken);
                autoPtr<entry> e = entry::New(is);
                dictionary d;
                d.add(e);
                return
                    autoPtr<Function1<scalar>>
                    (
                        Function1<scalar>::New(keyName, d)
                    );
            }
        };

        sigmaTable sigmas_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Dot product of face unit interface normal flux
        // and Face area vectors
        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Face unit interface normal flux
        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Function correcting contact angle
        void correctContactAngle
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            surfaceVectorField::Boundary& nHatb
        ) const;

        //- Simple expression for curvature
        tmp<volScalarField> K
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Surface tension force
        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U
        ) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
