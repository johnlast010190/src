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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::fv::interRegionHeatTransferModel

Group
    grpFvOptionsSources

Description
    Base class for inter region heat exchange. The derived classes must
    provide the heat transfer coeffisine (htc) which is used as follows
    in the energy equation.

     \f[
        -htc*Tmapped + Sp(htc, T)
     \f]

\*---------------------------------------------------------------------------*/

#ifndef interRegionHeatTransferModel_H
#define interRegionHeatTransferModel_H

#include "interRegionOption/interRegionOption.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class interRegionHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class interRegionHeatTransferModel
:
    public interRegionOption
{
protected:

    // Protected data

        //- Name of the model in the neighbour mesh
        word nbrModelName_;

        //- Pointer to neighbour interRegionHeatTransferModel
        interRegionHeatTransferModel* nbrModel_;

        //- First iteration
        bool firstIter_;

        //- Time index - used for updating htc
        label timeIndex_;

        //- Heat transfer coefficient [W/m2/k] times area/volume [1/m]
        volScalarField htc_;

        //- Flag to activate semi-implicit coupling
        bool semiImplicit_;

        //- Name of temperature field; default = "T"
        word TName_;

        //- Name of neighbour temperature field; default = "T"
        word TNbrName_;

        //- Reporting y/n
        bool reportOn_;

        //- total energy source
        scalar totalEnergy_;

        //- target energy source
        scalar targetEnergy_;

        //- Flag to ensure energy conservation on both sides
        bool conservative_;


    // Protected member functions

        //- Set the neighbour interRegionHeatTransferModel
        void setNbrModel();

        //- Inherit correct from interRegionOption
        using interRegionOption::correct;

        //- Correct to calculate the inter-region heat transfer coefficient
        void correct();

        //- Return the mapped temperature field
        tmp<volScalarField> getTmapped
        (
            const volScalarField& T
        ) const;

        //- Interpolate field with nbrModel specified
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const interRegionHeatTransferModel& nbrModel,
            const Field<Type>& field
        ) const;

        //- Interpolate field without nbrModel specified
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const Field<Type>& field
        ) const;

        //- Interpolate field with nbrModel specified
        template<class Type>
        void interpolate
        (
            const interRegionHeatTransferModel& nbrModel,
            const Field<Type>& field,
            Field<Type>& result
        ) const;

        //- Interpolate field without nbrModel specified
        template<class Type>
        void interpolate
        (
            const Field<Type>& field,
            Field<Type>& result
        ) const;


public:

    //- Runtime type information
    TypeName("interRegionHeatTransferModel");


    // Constructors

        //- Construct from dictionary
        interRegionHeatTransferModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~interRegionHeatTransferModel();


    // Member Functions

        // Access

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the mapToMap pointer
            inline const meshToMesh& meshInterp() const;

            //- Return the heat transfer coefficient
            inline const volScalarField& htc() const;

            //- Return const access to the neighbour model
            inline const interRegionHeatTransferModel& nbrModel() const;

            //- Return access to the neighbour model
            inline interRegionHeatTransferModel& nbrModel();

            //- Source term to energy equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Calculate heat transfer coefficient
            virtual void calculateHtc(scalar deltaT) = 0;

            //- Return total energy source
            const scalar& totalEnergy()
            {
                return totalEnergy_;
            }

            //- Return target energy source
            const scalar& targetEnergy()
            {
                return targetEnergy_;
            }

            //- Calculate target energy source
            virtual scalar getTargetEnergy();


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

            //- Do reporting y/n
            inline bool reportOn()
            {
                return reportOn_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interRegionHeatTransferModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interRegionHeatTransferModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
