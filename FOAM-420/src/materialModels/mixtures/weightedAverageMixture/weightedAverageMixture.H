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
    (c) 2021-2023 Esi Ltd.

Class
    Foam::weightedAverageMixture

Description
    Homogeneus mixture model.

SourceFiles
    weightedAverageMixture.C

\*---------------------------------------------------------------------------*/

#ifndef weightedAverageMixture_H
#define weightedAverageMixture_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "mixtures/phaseVolumeFractions/phaseVolumeFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class weightedAverageMixture Declaration
\*---------------------------------------------------------------------------*/

class weightedAverageMixture
:
    public materialModel,
    public weightedAverageMixtureModel,
    public vectorWeightedAverageMixtureModel,
    public tensorWeightedAverageMixtureModel
{
    // Private member data

        //- Reference to the class providing mixing fractions
        const volumeMassFractions& frac_;

        //- Lookup or construct fractions class
        const volumeMassFractions& lookupOrConstructBase();

public:


    // Constructors

        //- Construct from dictionary
        weightedAverageMixture
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& name,
            const word& phaseName,
            const word& specieName
        );

        //- New constructor
        static autoPtr<weightedAverageMixture> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );



    //- Destructor
    virtual ~weightedAverageMixture()
    {}

    // Member functions

        //- Return the instantiated type name
        TypeName("weightedAverageMixture");


    // Member functions

        //- Internal field mixing
        template<class Type>
        tmp<Field<Type>> weightedAverageMixtureTypeInternal
        (
            const UPtrList<baseModels<Type>>& mod
        ) const;

        //- Geometric field mixing
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        weightedAverageMixtureTypeGeometric
        (
            const UPtrList<baseModels<Type>>& mod
        ) const;

        //- Patch mixing
        template<class Type>
        tmp<Field<Type>> weightedAverageMixtureTypePatch
        (
            const UPtrList<baseModels<Type>>& mod,
            const label patchi
        ) const;

        //- Cell mixing
        template<class Type>
        Type weightedAverageMixtureTypeCell
        (
            const UPtrList<baseModels<Type>>& mod,
            const label celli
        ) const;

        bool isweightedAverageMixtureConst() const override
        {
            // Ignore the edge case in which all constituents are constant and
            // equal, as this is not common enough to warrant the increase in
            // efficiency gained by treating it as constant
            return false;
        }

        bool isvectorWeightedAverageMixtureConst() const override
        {
            return false;
        }

        bool istensorWeightedAverageMixtureConst() const override
        {
            return false;
        }

        declareScalarFieldFunctions(weightedAverageMixture);

        declareVectorFieldFunctions(vectorWeightedAverageMixture);

        declareTensorFieldFunctions(tensorWeightedAverageMixture);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Cast to scalar model
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Cast to vector model
        virtual baseModels<vector>* castVectorModel
        (
            const word& modelName
        ) override;

        //- Cast to tensor model
        virtual baseModels<tensor>* castTensorModel
        (
            const word& modelName
        ) override;

        //- Cast to tensor model
        virtual bool mixture() const override
        {
            return true;
        }

        //- Re-read data
        virtual bool read() override
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "weightedAverageMixtureTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
