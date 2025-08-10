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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2010-2017, 2020 Esi Ltd

Class
    Foam::radiation::fvDOM

Group
    grpRadiationModels

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter.

    Available absorption models:
    - constantAbsorptionEmission
    - greyMeanAbsoprtionEmission
    - wideBandAbsorptionEmission

Usage
    \verbatim
        fvDOMCoeffs
        {
            nPhi        4;          // azimuthal angles in PI/2 on X-Y.
                                    //(from Y to X)
            nTheta      0;          // polar angles in PI (from Z to X-Y plane)
            convergence 1e-3;       // convergence criteria for radiation
                                    //iteration
            maxIter     4;          // maximum number of iterations
             meshOrientation    (1 1 1); //Mesh ortientation used for 2D and 1D
        }

        solverFreq   1; // Number of flow iterations per radiation iteration
    \endverbatim

    The total number of solid angles is  4*nPhi*nTheta in 3-D.

    Operating modes:
    - 1-D:
      - ray directions are on X, Y or Z
      - \c nPhi and \c nTheta entries are ignored
      - \c meshOrientation vector can be used for any other 1-D direction.
    - 2-D:
      - ray directions are on X-Y, X-Z or Y-Z planes
      - only the \c nPhi entry is considered
      - \c meshOrientation vector can be used for non-aligned planes
        specifying the plane normal vector.
    - 3-D:
      - rays geberated in 3-D using the \c nPhi and \c nTheta entries
      - \c meshOrientation vector is not applicable.

Update:
    nRayInput switch provides an alternative method for angular discretisation.

    A total number of solid angle is specified with nPhi and nTheta calculated
    following empirical rules. nPhi and nTheta are chosen to give a total
    number of solid angles as close as possible to the value specified by
    the user.

        fvDOMCoeffs
        {
            nRayInput true;
            nRay 32;
        }

    The total number of solid angles is nPhi*nTheta in 3-D and nPhi in 2D.
    Note that theta [0, PI] and phi [0, 2*PI].

    This radiation model also supports solar/environmental radiation in the
    form of solarLoad and domSolar. See
    radiation/radiationModels/solarLoad/solarLoad.H
    radiation/radiationModels/domSolar/domSolar.H
SourceFiles
    fvDOM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelfvDOM_H
#define radiationModelfvDOM_H

#include "radiationModels/fvDOM/radiativeIntensityRay/radiativeIntensityRay.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "fvMatrices/fvMatrices.H"
#include "cfdTools/general/explicitScalarWaveSolve/explicitScalarWaveSolve.H"
#include "radiationModels/solarLoad/solarLoad.H"
#include "radiationModels/domSolar/domSolar.H"
#include "submodels/boundaryRadiationProperties/boundaryRadiationProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class fvDOM Declaration
\*---------------------------------------------------------------------------*/

class fvDOM
:
    public radiationModel
{

public:

        enum solarLoadType
        {
            sltViewFactor,
            sltDOM,
            sltnone
        };

        static const NamedEnum<solarLoadType, 3> solarLoadTypeNames_;

private:
    // Private data


        //- Incident radiation  [W/m2]
        volScalarField G_;

        //- Total (net) radiative heat flux [W/m2]
        volScalarField qr_;

        //- Incident (absorbed) radiative heat flux [W/m2]
        volScalarField qin_;

        //- Emmited radiative heat flux [W/m2]
        // volScalarField qem_;

        //- Irradiation heat flux [W/m2]
        volScalarField qg_;

        //- Total absorption coefficient [1/m]
        volScalarField a_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        // Calculate nTheta_ and nPhi_ based on
        // total number of rays provided as input
        Switch nRayInput_;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nLambda_;

        //- Wavelength total absorption coefficient [1/m]
        PtrList<volScalarField> aLambda_;

        //- Black body
        blackBodyEmission blackBody_;

        //- List of pointers to radiative intensity rays
        PtrList<radiativeIntensityRay> IRay_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        label maxIter_;

        //- Initial convergence criterion
        scalar initConvergence_;

        //- Initial maximum number of iterations
        scalar initMaxIter_;

        //- List of cached wave solver coefficients for rays
        List<PtrList<waveSolverCoeffs>>rayWaveCoeffs_;

        //- Maximum omega weight
        scalar omegaMax_;

        //- Solar Load model type selector
        solarLoadType solarLoadMode_;

        //- View factor solar load radiation model
        autoPtr<solarLoad> vfSolarLoad_;

        //- fvDOM solar load radiation model
        autoPtr<domSolar> domSolarLoad_;


        //- Mesh orientation vector
        vector meshOrientation_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        fvDOM(const fvDOM&);

        //- Disallow default bitwise assignment
        void operator=(const fvDOM&);

        //- Update nlack body emission
        void updateBlackBodyEmission();


public:

    //- Runtime type information
    TypeName("fvDOM");


    // Constructors

        //- Construct from components
        fvDOM(const volScalarField& T);

        //- Construct from components
        fvDOM(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~fvDOM();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayIdLambdaId
            (
                const word& name,
                label& rayId,
                label& lambdaId
            ) const;

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;

        //- calculate patch heat flux due to immitted radiation
            tmp<scalarField> emittedRadiantIntensity
            (
                label patchI,
                const scalarField& Tp
            ) const;

        // Access

            //- Ray intensity for rayI
            inline const radiativeIntensityRay& IRay(const label rayI) const;

            //- Ray intensity for rayI and lambda bandwidth
            inline const volScalarField& IRayLambda
            (
                const label rayI,
                const label lambdaI
            ) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of wavelengths
            inline label nLambda() const;

            //- Const access to total absorption coefficient
            inline const volScalarField& a() const;

            //- Const access to wavelength total absorption coefficient
            inline const volScalarField& aLambda(const label lambdaI) const;

            //- Const access to incident radiation field
            inline const volScalarField& G() const;

            //- Const access to total (net) radiative heat flux field
            inline const volScalarField& qr() const;

            //- Const access to incident (absorbed) radiative heat flux field
            inline const volScalarField& qin() const;

            //- Const access to emitted radiative heat flux field
            //inline const volScalarField& qem() const;

            //- Const access to irradiation heat flux field
            inline const volScalarField& qg() const;

            //- Const access to black body
            inline const blackBodyEmission& blackBody() const;

            //- Const access to wave solver coefficients
            inline const waveSolverCoeffs& rayWaveCoeffs
            (
                const label lambdaI,
                const label rayId
            ) const;

            //- Return omegaMax
            inline scalar omegaMax() const;

            //- Return meshOrientation
            inline vector meshOrientation() const;

            inline const boundaryRadiationProperties&
            boundaryProperties() const;

            const domSolar& getDomSolarObj() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radiationModels/fvDOM/fvDOM/fvDOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
