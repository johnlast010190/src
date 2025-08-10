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
    (c) ESI Ltd, 2019.

Description :
  This file declares the objects and functions for the barotropic cavitation model
  implementation. FoamVOF is a general-purpose multi-fluid module developed
  under the framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/
#ifndef barocav_H
#define barocav_H

namespace Foam
{
    class BaroCav : public FoamVOF
    {
        public:
            scalar CoNum;
            scalar meanCoNum ;
            scalar acousticCoNum;
            scalar maxAcousticCo;
            bool correctPhi_;
            dimensionedScalar psil;
            dimensionedScalar rholSat;
            dimensionedScalar psiv;
            dimensionedScalar pSat;
            dimensionedScalar rhovSat;
            dimensionedScalar rhol0;
            dimensionedScalar rhoMin;

            bool correctPhi()
            {
                return correctPhi_;
            }
            void correctPhi(pimpleControl &pimple);

            autoPtr<barotropicCompressibilityModel> psiModel;

            BaroCav(const Time& time, const dynamicFvMesh& mesh,const MatProp &prop);

            void CourantNo();
            void createFields(pimpleControl &pimple);

            void rhoEqnSolve();
            void correctAlphav();

            autoPtr<wordList> pcorrTypes;


            void updateFields
            (
                pimpleControl &pimple,
                incompressible::turbulenceModel &turbulence,
                fv::options& fvOptions
            );


            void UEqnSolve
            (
                fvVectorMatrix& UEqn,
                pimpleControl &pimple,
                incompressible::turbulenceModel& turbulence
            );

            void pEqnSolve
            (
                fvVectorMatrix& UEqn,
                pimpleControl &pimple,
                incompressible::turbulenceModel& turbulence
            );


            const volScalarField& psi()
            {
                return psiModel().psi();
            }
            //alpha1 is vapour

            volScalarField& alpha1();
            volScalarField& alpha2();

            volScalarField& alphav()
            {
                return mixture2().alpha1();
            }

            volScalarField& alphal()
            {
                return mixture2().alpha2();
            }

    };

}// End namespace Foam

#endif
