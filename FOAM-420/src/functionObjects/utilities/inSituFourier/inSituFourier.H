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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2015 Esi Ltd.

Class
     Foam::functionObjects::inSituFourier

Group
    grpUtilitiesFunctionObjects

Description

Usage
    \verbatim
    inSituFourier
    {
        type                inSituFourier;
        libs                ("libutilityFunctionObjects.so");

        writeControl        writeTime;

        //Field to use for Fourier; defaults to p
        field           p;

        //Compute on internal field (default = yes)
        doVolume        yes;

        //Only compute on specified patches (optional)
        //patches         ("fixedWalls");

        //Starting time of Fourier analysis (optional)
        //(default = start of initial analysis)
        timeStart             0.0;

        //Reference value for dB calculation (default as shown)
        dBRefVal        dBRefVal [ 1 -1 -2 0 0 0 0 ] 20e-6;

        //Automatically generate frequency bands
        generateFrequencyBands   yes;

        //Number of frequency bands to generate per octave (default = 3)
        bandsPerOctave  12;

        //Lowest centre frequency for auto-generated bands
        fl              315;

        //Maximum centre frequency for auto-generated bands
        fU              16000;

        //Manual list of frequency bands (if generateFrequencyBands = no)
        frequencyBands ((100  200)(200 400));

        //Centre frequencies (if generateFrequencyBands = no) (optional).
        //Defaults to geometric mean of band min/max
        centreFrequencies (150 300);
    }
    \endverbatim

SourceFiles
    inSituFourier.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_inSituFourier_H
#define functionObjects_inSituFourier_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class inSituFourier Declaration
\*---------------------------------------------------------------------------*/

class inSituFourier
:
    public fvMeshFunctionObject
{
    // Private data

        autoPtr<dimensionedScalar> PhiRescaling_;

        // Read from dictonary

            //- input
            word fieldName_;

            //- start time of Fourier analysis
            scalar t0_;

            //- calculate on internal field
            Switch doVolume_;

            //- Optional calculation on selected patches
            labelHashSet patchSet_;

            //- numerical constants
            autoPtr<dimensionedScalar> dBRefVal_; //ref field value for dB calc
            scalar ecf_; //energy conversion factor

            //- frequency bands
            List<scalar> flow_;
            List<scalar> fm_;
            List<scalar> fhigh_;

        //- Computed values

            PtrList<volScalarField> avgPower_;
            PtrList<volScalarField> dB_;
            PtrList<volScalarField> integralPhiSinWinDt_;
            PtrList<volScalarField> integralPhiCosWinDt_;
            PtrList<uniformDimensionedScalarField> integralSinWinDt_;
            PtrList<uniformDimensionedScalarField> integralCosWinDt_;
            PtrList<volScalarField> avgPhi_;


    // Private Member Functions

        void createFields();

        void deactivateFields(volScalarField& f);

        void doIntegration
        (
            const label bandNo,
            const scalarField& Phi,
            scalarField& integralPhiSinWinDt,
            scalarField& integralPhiCosWinDt,
            scalar& integralSinWinDt,
            scalar& integralCosWinDt,
            scalarField& avgPhi,
            scalarField& avgPower,
            scalarField& dB
        );

        //- Disallow default bitwise copy construct and assignment
        inSituFourier(const inSituFourier&);
        void operator=(const inSituFourier&);


public:

    //- Runtime type information
    TypeName("inSituFourier");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        inSituFourier
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    // Destructor

        virtual ~inSituFourier();


    // Member Functions

        //- Update solver
        void calculate();

        //- Read the parameters
        virtual bool read(const dictionary&);

        //- Execute fourier calculation
        virtual bool execute();

        //- Sample and write
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
