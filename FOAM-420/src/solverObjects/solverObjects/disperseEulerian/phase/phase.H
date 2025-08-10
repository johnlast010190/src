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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::decoupledEulerian::phase

Description

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"

#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

namespace Foam
{

namespace decoupledEulerian
{

/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class phase
:
    public regIOobject
{
protected:

    // Protected data
        const word name_;
        const dictionary dict_;

        //- Material data continuous phase
        const dimensionedScalar alphac_;
        const objectRegistry& obr_;
        const fvMesh& mesh_;

        mutable autoPtr<volScalarField> rhoc_;
        mutable autoPtr<volScalarField> muc_;
        dimensionedScalar rhocVal_;
        dimensionedScalar mucVal_;
        dimensionedScalar kappac_;   // heat transfer
        scalar Pr_;                  // heat transfer
        dimensionedScalar Tsat_;     // phase change
        dimensionedScalar L_;        // phase change
        dimensionedScalar Dc_;       // physisorption

        //- Material data disperse phase
        autoPtr<volScalarField> alphad_;
        autoPtr<volVectorField> Ud_;
        dimensionedScalar rhod_;
        dimensionedScalar diam_;
        dimensionedScalar Td_;       // heat transfer
        dimensionedScalar Md_;       // phase change
        dimensionedScalar dHevap_;   // phase change


    // Protected member functions
        void getRhod();
        void getDiam();
        void updateRhoc() const;
        void updateMuc() const;


public:

    //- Runtime type information
    TypeName("phase");


    // Constructors

        // Construct from dictionary
        phase
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr,
            const fvMesh& mesh,
            const bool registerObject = true
        );


    //- Destructor
    virtual ~phase();


    // Access Functions

        //- Phase name
        const word& name() const {return name_;}

        //- Phase dict
        const dictionary& dict() const {return dict_;}

        //- Continuous phase density
        const volScalarField& rhoc() const
        {
            // update rhoc from turbulence model
            updateRhoc();
            return rhoc_();
        }

        //- Continuous phase viscosity
        const volScalarField& muc() const
        {
            // update muc from turbulence model
            updateMuc();
            return muc_();
        }

        //- Continuous phase thermal conductivity
        const dimensionedScalar& kappac() const {return kappac_;}

        //- Continuous phase Prandtl number
        const scalar& Pr() const {return Pr_;}

        //- Continuous phase saturation temperature
        const dimensionedScalar& Tsat() const {return Tsat_;}

        //- Latent heat
        const dimensionedScalar& L() const {return L_;}

        //- Continuous phase diffusion coefficient
        const dimensionedScalar& Dc() const {return Dc_;}

        //- Disperse phase density
        const dimensionedScalar& rhod() const {return rhod_;}

        //- Disperse phase diameter
        const dimensionedScalar& diam() const {return diam_;}

        //- Disperse phase temperature
        const dimensionedScalar& Td() const {return Td_;}

        //- Disperse phase Molar mass
        const dimensionedScalar& Md() const {return Md_;}

        //- Disperse phase evaporation enthalpy
        const dimensionedScalar& dHevap() const {return dHevap_;}

        //- Disperse phase Reynolds number
        tmp<volScalarField> Re() const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr() const;

        //- Relative velocity
        tmp<volVectorField> Ur() const;

        //- Access phase volume fraction
        const volScalarField& alphad() const {return alphad_();}
        volScalarField& alphad() {return alphad_();}

        //- Access continuous phase velocity
        tmp<volVectorField> Uc() const
        {
            const turbulenceModel& turb =
                obr_.lookupObject<turbulenceModel>(turbulenceModel::propertiesName);
            return turb.U();
        }

        //- Access disperse phase velocity
        const volVectorField& Ud() const {return Ud_();}
        volVectorField& Ud() {return Ud_();}

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;

        //- Read data from dictionary
        void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decoupledEulerian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
