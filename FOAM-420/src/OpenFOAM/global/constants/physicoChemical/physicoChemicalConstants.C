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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2014 blueCAPE Lda

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.



    Modifications made:
      - These changes are technically from Symscape's own patches for Windows,
        circa 2011.
      - The adaptation derived from the patches for blueCFD 2.1, adjusted to
        OpenFOAM 2.2.


\*---------------------------------------------------------------------------*/

#include "global/constants/mathematical/mathematicalConstants.H"
#include "global/constants/universal/universalConstants.H"
#include "global/constants/electromagnetic/electromagneticConstants.H"
#include "global/constants/physicoChemical/physicoChemicalConstants.H"

#include "global/constants/dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{

const char* const physicoChemical::group = "physicoChemical";

defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::R,
    dimensionedScalar
    (
        "R",
        physicoChemical::NA*physicoChemical::k
    ),
    constantphysicoChemicalR,
    "R"
);


defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::F,
    dimensionedScalar
    (
        "F",
        physicoChemical::NA*electromagnetic::e
    ),
    constantphysicoChemicalF,
    "F"
);


#if defined( WIN32 ) || defined( WIN64 )
// Note: cannot use dimless etc. since not guaranteed to be constructed
defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::sigma,
    dimensionedScalar
    (
        "sigma",
        (
            Foam::dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),    //Foam::dimless,
                Foam::sqr(mathematical::pi)/60.0
            )
        *Foam::pow4(physicoChemical::k)
        /(pow3(universal::hr)*sqr(universal::c))
       )
#ifdef FOAM_SP
       .dimensions(),
       // Assuming this is the Stefan-Boltzmann constant
       // http://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law
       // Single precision can't handle the pow4(k),
       // where k is Boltzmann constant = 1.3806488e-23
       5.6704e-8f
#endif
    ),
    constantphysicoChemicalsigma,
    "sigma"
);

#else
defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::sigma,
    // Avoid overflow in single precision
    dimensionedScalar
    (
        "sigma",
        Foam::dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),    //Foam::dimless,
            Foam::sqr(mathematical::pi)/60.0
        )
       *sqr(physicoChemical::k/(universal::hr*universal::c))
       *physicoChemical::k/universal::hr*physicoChemical::k
    ),
    constantphysicoChemicalsigma,
    "sigma"
);
#endif

defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::b,
    dimensionedScalar
    (
        "b",
        (universal::h*universal::c/physicoChemical::k)
       /Foam::dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),    //Foam::dimless
            4.965114231
        )
    ),
    constantphysicoChemicalb,
    "b"
);


defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::c1,
    dimensionedScalar
    (
        "c1",
        Foam::dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),    //Foam::dimless,
            mathematical::twoPi
        )
       *universal::h*Foam::sqr(universal::c)
    ),
    constantphysicoChemicalc1,
    "c1"
);


defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::c2,
    dimensionedScalar
    (
        "c2",
        universal::h*universal::c/physicoChemical::k
    ),
    constantphysicoChemicalc2,
    "c2"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
