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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::groupOfHumansHeatFluxFvPatchScalarField

Description
    Empirical mass flux boundary for heat rejection.

SourceFiles
    groupOfHumansHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef groupOfHumansHeatFluxFvPatchScalarFields_H
#define groupOfHumansHeatFluxFvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "turbulentTransportModels/derivedFvPatchFields/turbulentHeatFluxTemperature/turbulentHeatFluxTemperatureFvPatchScalarField.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
     Class groupOfHumansHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class groupOfHumansHeatFluxFvPatchScalarField
:
    public turbulentHeatFluxTemperatureFvPatchScalarField
{
public:

    // group enumeration
    enum eGroup
    {
        mA,
        fA,
        mC,
        fC
    };

    // human struct
    struct demographic
    {
        public:
        scalar N;
        scalar height;
        scalar weight;
        scalar age;

        demographic(const dictionary& dict)
        :
            N(readScalar(dict.lookup("number"))),
            height(readScalar(dict.lookup("height"))),
            weight(readScalar(dict.lookup("weight"))),
            age(dict.lookupOrDefault<scalar>("age", 0.0))
        {}

        demographic()
        :
            N(0.0),
            height(0.0),
            weight(0.0),
            age(0.0)
        {}

        void addTo(dictionary& dict) const
        {
            dict.add("number", N);
            dict.add("height", height);
            dict.add("weight", weight);
            dict.add("age", age);
        }

    };

protected:

    static const NamedEnum<eGroup, 4> eGroupNames_;


private:

    // Input data

        //- population
        List<demographic> people_;

        //- Standard metabolic rate
        scalar Mstd_;

        //max distance from surface for local cells
        scalar dlocal_;

        // mean far wall emissivity
        scalar emissivityFarWall_;

        // mean local emissivity
        scalar emissivityLocal_;

        // StephanBoltzmann constant
        scalar sigmaSB_;

    // Derived data
        //- populatiuon stats
        //- Effective weight
        scalar Weff_;

        //- Effective height
        scalar Heff_;

        //- Effective area
        scalar Aeff_;

        //- Actual area
        scalar Areal_;

        //- cellSet for calculating local mean temperature and velocity
        cellSet nearCells_;

        //- total volume of near cells
        scalar nearVolume_;

        //- total number of people
        scalar Ntotal_;

    // Demand driven data

        //- Effective metabolic rate
        mutable autoPtr<scalar> MeffPtr_;

        //- average local pressure
        mutable autoPtr<scalar> PlocalPtr_;

        //- average local temperature
        mutable autoPtr<scalar> TlocalPtr_;

        //- effective room temperature
        mutable autoPtr<scalar> TroomPtr_;


    // Private member functions

        //- collect all cells within specified distance d in cellSet nearCells_
        void findLocalCells(const fvMesh& mesh, scalar d);

        //- surface area of human calculated from height and weight
        scalar humanSurfaceArea(scalar W,scalar H) const;

        //- calculate effective height, weight and area for population
        void staticStatisticsCalculator();

        //- metabolic heat flux [W/m2]
        scalar humanHeatFlux() const;

        //- effective metabolic rate for a group of humans
        scalar effectiveMetabolicRate() const;

        //- mean radiant temperature
        scalar meanRadiantTemperature() const;

        //- makePeople
        List<demographic> makePeople(const dictionary&) const;

        //- calculate effective room temperature
        scalar effectiveRoomtemperature() const;

        //- mean value of field magnitude in cellSet nearCells_
        template<class volField>
        scalar calcLocalMean(word f) const;


public:

    //- Runtime type information
    TypeName("incompressible::groupOfHumansHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        groupOfHumansHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        groupOfHumansHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  groupOfHumansHeatFluxFvPatchScalarField onto
        //  a new patch
        groupOfHumansHeatFluxFvPatchScalarField
        (
            const groupOfHumansHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        groupOfHumansHeatFluxFvPatchScalarField
        (
            const groupOfHumansHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new groupOfHumansHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        groupOfHumansHeatFluxFvPatchScalarField
        (
            const groupOfHumansHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new groupOfHumansHeatFluxFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // Access

            //- info needed by dependent physiological boundaries

            scalar Plocal() const;

            scalar Tlocal() const;

            scalar Troom() const;

            scalar Meff() const;

            scalar Aeff() const
            {
                return Aeff_;
            }

            scalar Ntotal() const
            {
                return Ntotal_;
            }

            scalar Areal() const
            {
                return Areal_;
            }


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
