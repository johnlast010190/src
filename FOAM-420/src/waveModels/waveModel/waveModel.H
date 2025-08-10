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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2015 IH-Cantabria

Class
    Foam::waveModel

Description
    Base class for waveModels

\*---------------------------------------------------------------------------*/

#ifndef waveModel_H
#define waveModel_H

#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/dictionary/dictionary.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "fields/Fields/vectorField/vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                          Class waveModel Declaration
\*---------------------------------------------------------------------------*/

class waveModel
:
    public refCount,
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the polyPatch
        const polyPatch& patch_;

        //- Gravity
        const vector& g_;

        //- Name of velocity field, default = "U"
        word UName_;

        //- Name of phase fraction field, default = "alpha"
        word alphaName_;

        //- Rotation tensor from global to local system
        tensor Rgl_;

        //- Rotation tensor from local to global system
        tensor Rlg_;

        //- Number of paddles
        label nPaddle_;

        //- Paddle x co-ordinates / [m]
        scalarField xPaddle_;

        //- Paddle y co-ordinates / [m]
        scalarField yPaddle_;

        //- Addressing from patch face index to paddle index
        labelList faceToPaddle_;

        //- Patch face centre z co-ordinates / [m]
        scalarField z_;

        //- Overall (point) span in z-direction / [m]
        scalar zSpan_;

        //- Minimum z (point) height per patch face / [m]
        scalarField zMin_;

        //- Maximum z (point) height per patch face / [m]
        scalarField zMax_;

        //- Minimum z reference level / [m]
        scalar zMin0_;

        //- Reference water depth / [m]
        scalar waterDepthRef_;

        //- Initial depth / [m]
        scalar initialDepth_;

        //- Time index used for updating
        label currTimeIndex_;

        //- Active wave absorption switch
        bool activeAbsorption_;


        // Current values

            //- Velocity field
            vectorField U_;

            //- Wave indicator field
            scalarField alpha_;


    // Protected Member Functions

        //- Initialise
        virtual void initialiseGeometry();

        //- Water level
        virtual tmp<scalarField> waterLevel() const;

        //- Return the time scaling coefficient
        virtual scalar timeCoeff(const scalar t) const = 0;

        //- Set the water level
        virtual void setLevel
        (
            const scalar t,
            const scalar tCoeff,
            scalarField& level
        ) const = 0;

        //- Calculate the wave model velocity
        virtual void setVelocity
        (
            const scalar t,
            const scalar tCoeff,
            const scalarField& level
        ) = 0;

        //- Set the alpha field based on the water level
        virtual void setAlpha(const scalarField& level);

        //- Set the paddle coverage fraction and reference height
        virtual void setPaddlePropeties
        (
            const scalarField& level,
            const label facei,
            scalar& fraction,
            scalar& z
        ) const;


public:

    //- Runtime type information
    TypeName("waveModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            waveModel,
            patch,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const polyPatch& patch
            ),
            (dict, mesh, patch)
        );


    //- Constructor
    waveModel
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );


    // Selectors

        //- Return a reference to the selected wave model
        static autoPtr<waveModel> New
        (
            const word& dictName,
            const fvMesh& mesh,
            const polyPatch& patch
        );

        //- Lookup waveModel from database, or create new
        static tmp<waveModel> lookupOrCreate
        (
            const polyPatch& patch,
            const fvMesh& mesh,
            const word& waveDictName
        );


    //- Destructor
    virtual ~waveModel();

    //- Dictionary name
    static const word dictName;


    // Public Member Functions

        //- Utility function to construct the model name
        static word modelName(const word& patchName);

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Return the latest wave velocity prediction
        virtual const vectorField& U() const;

        //- Return the latest wave indicator field prediction
        virtual const scalarField& alpha() const;

        //- Correct the model for time, t[s]
        virtual void correct(const scalar t);

        //- Info
        virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
