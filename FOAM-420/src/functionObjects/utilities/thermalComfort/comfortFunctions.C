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
  (c) 2009, Peter Broede
  (c) 2010-2011, Esi Ltd

\*---------------------------------------------------------------------------*/

#include "thermalComfort/comfortFunctions.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(comfortFunctions, 0);
}

Foam::scalar Foam::comfortFunctions::wetBulb
(
    const scalar& temperature,
    const scalar& relativeHumidity
 )
{
        //Convert to Celsius
    scalar Ta = temperature - 273.15;
    scalar RH = relativeHumidity;

    scalar wb = Ta*atan(scalar(0.151977)
                        *pow(RH+scalar(8.313659),0.5))
        + atan(Ta+RH)
        - atan(RH - scalar(1.676331))
        + scalar(0.00391838)*pow(RH,1.5)
        * atan(0.023101*RH) - 4.686035;

    return wb + 273.15;
}


    // Calculate water vapour pressure (Torr)
Foam::scalar Foam::comfortFunctions::svpTorr
(
    const scalar Ta
 )
{
    return 0.007500637554*svp(Ta);
}


    // Calculate water vapour pressure (Pa)
Foam::scalar Foam::comfortFunctions::svp
(
    const scalar Ta
 )
{
    List<scalar> g(8);

    g[0] = -2.8365744E+03;
    g[1] = -6.028076559E+03;
    g[2] = 1.954263612E+01;
    g[3] = -2.737830188E-02;
    g[4] = 1.6261698E-05;
    g[5] = 7.0229056E-10;
    g[6] = -1.8680009E-13;
    g[7] = 2.7150305;

    scalar e = g[7]*log(Ta);

    for (int i = 0; i < 7; i++)
    {
        e += g[i]*pow(Ta, i-2);
    }
    e = exp(e);

    return e;
}


Foam::scalar Foam::comfortFunctions::TOp
(
    const scalar& va,
    const scalar& temperature,
    const scalar& radiantTemperature
 )
{
    scalar Ta = temperature - 273.15;
    scalar Tmrt = radiantTemperature - 273.15;

    scalar A = 0.7;
    if (va < 0.2)
    {
        A = 0.5;
    }
    else if (va < 0.6)
    {
        A = 0.6;
    }

    return
        (
            A*Ta + (1.0 - A) * Tmrt
         );
}


Foam::scalar Foam::comfortFunctions::PPD
(
    const scalar& PMV
 )
{
    return
        (
            100.0 - 95 * Foam::exp
            (
                -0.03353*pow(PMV,4) - 0.2179 * pow(PMV,2)
             )
         );
}


Foam::scalar Foam::comfortFunctions::PMV
(
    const scalar& va,
    const scalar& temperature,
    const scalar& radiantTemperature,
    const scalar& relativeHumidity
 )
{
    scalar P1(0), P2(0), P3(0), P4(0), P5(0);
    scalar XN(0), XF(0), HCN(0), HCF(0), HC(0), FCL(0), EPS(0);
    scalar HL1(0), HL2(0), HL3(0), HL4(0), HL5(0), HL6(0);
    scalar TCL(0), TS(0), TCLA(0);
    label N;

        //Water vapour pressure (Pa)
    scalar wvp = svp(temperature)*relativeHumidity/100.0;

    if (clo_ <= 0.078) // ok
    {
        FCL = 1 + 1.29 * clo_;
    } // ok
    else
    {
        FCL = 1.05 + 0.645 * clo_; // clothing area factor
    }

    HCF = 12.1 * Foam::sqrt(va);
    TCLA = temperature + (35.5 - (temperature - 273.15))
        / (3.5 * ( 6.45*clo_ + 0.1));

    XN = TCLA / 100.0; // ok

    P1 = clo_ * FCL; // ok
    P2 = P1 * 3.96; // ok
    P3 = P1 * 100.0; // ok
    P4 = P1 * temperature; // ok
    P5 = 308.7 - 0.028 * (met_ - work_)
        + P2 * Foam::pow(radiantTemperature/100.0, 4);
    XF = XN;
    EPS = 0.00015;

    N=0;
    do
    {
        N++;
        XF = (XF + XN)/2; // stop criteria by iteration
        HCN = 2.38 * Foam::pow( mag( 100.0 * XF - temperature ), 0.25);

        if (HCF>HCN) { HC = HCF; } else { HC = HCN; }; // ok
        XN = (P5 + P4 * HC - P2 * Foam::pow(XF,4) )
            / (100.0 + P3 * HC);
        if (N>150) { break; };
    } while (mag(XN-XF)>EPS);

    TCL = 100.0 * XN - 273.15; // surface temperature of clothing
    HL1 = 3.05 * 0.001
        * (5733 - 6.99*(met_ - work_) - wvp);
        // heat loss diff. through skin
    if ((met_ - work_) > 58.15)
    {
        HL2 = 0.42 * ( (met_ - work_) - 58.15);
    }
    else
    {
        HL2 = 0;
    } // heat loss by sweating (comfort)

        // latent respiration heat loss
    HL3 = 1.7 * 0.00001 * met_ * (5867 - wvp);
        // dry respiration heat loss
    HL4 = 0.0014 * met_ * (34 - (temperature - 273.15) );
        // heat lose by radiation
    HL5 = 3.96 * FCL * (Foam::pow(XN,4)
                        - Foam::pow( (radiantTemperature/100.0),4)) ;
        // heat lose by convection
    HL6 = FCL * HC * (TCL - (temperature - 273.15));
        // thermal sensation trans coeff
    TS = 0.303 * Foam::exp(-0.036 * met_ ) + 0.028;

        // Berechne PMV-Wert
    return
        (
            TS * ( (met_ - work_)
                   - HL1 - HL2 - HL3 - HL4 - HL5 - HL6 )
         );
}


Foam::scalar Foam::comfortFunctions::DR
(
    const scalar& va,
    const scalar& temperature,
    const scalar& Tu_
)
{
    scalar DR;

    if (va >= 0.05)
    {
        DR = ( 34 - (temperature - 273.15) ) *
            ( Foam::pow( va - 0.05 ,0.62 )
              * ( (0.37 * va * Tu_) + 3.14) );
    }
    else
    {
        DR = (34 - temperature - 273.15) *
            Foam::pow( 0 ,0.6223) * ((0.37 * 0.05 * Tu_) + 3.14);
    }
    if (DR > 100)
    {
        DR = 100;
    }

    return DR;
}

Foam::scalar Foam::comfortFunctions::TU
(
    const scalar& va,
    const scalar& k,
    const scalar upperLimitTu,
    const scalar lowerLimitTu
)
{
    return max((min((100.0*sqrt(k*2/3)/va),upperLimitTu)),lowerLimitTu);
}


Foam::scalar Foam::comfortFunctions::UTCI
(
    const scalar& va,
    const scalar& temperature,
    const scalar& radiantTemperature,
    const scalar& relativeHumidity
 )
{
        //Convert to Celsius
    scalar Ta = temperature - 273.15;
    scalar Tmrt = radiantTemperature - 273.15;

        //Water vapour pressure (Pa)
    scalar wvp = svp(temperature)*relativeHumidity/100.0;

    scalar DTmrt = Tmrt-Ta;
        //Convert vapour pressure from Pa to kPa
    scalar kPa = wvp / 1000;

    return Ta +
        ( 6.07562052E-01 )   +
        ( -2.27712343E-02 ) * Ta +
        ( 8.06470249E-04 ) * Ta*Ta +
        ( -1.54271372E-04 ) * Ta*Ta*Ta +
        ( -3.24651735E-06 ) * Ta*Ta*Ta*Ta +
        ( 7.32602852E-08 ) * Ta*Ta*Ta*Ta*Ta +
        ( 1.35959073E-09 ) * Ta*Ta*Ta*Ta*Ta*Ta +
        ( -2.25836520E+00 ) * va +
        ( 8.80326035E-02 ) * Ta*va +
        ( 2.16844454E-03 ) * Ta*Ta*va +
        ( -1.53347087E-05 ) * Ta*Ta*Ta*va +
        ( -5.72983704E-07 ) * Ta*Ta*Ta*Ta*va +
        ( -2.55090145E-09 ) * Ta*Ta*Ta*Ta*Ta*va +
        ( -7.51269505E-01 ) * va*va +
        ( -4.08350271E-03 ) * Ta*va*va +
        ( -5.21670675E-05 ) * Ta*Ta*va*va +
        ( 1.94544667E-06 ) * Ta*Ta*Ta*va*va +
        ( 1.14099531E-08 ) * Ta*Ta*Ta*Ta*va*va +
        ( 1.58137256E-01 ) * va*va*va +
        ( -6.57263143E-05 ) * Ta*va*va*va +
        ( 2.22697524E-07 ) * Ta*Ta*va*va*va +
        ( -4.16117031E-08 ) * Ta*Ta*Ta*va*va*va +
        ( -1.27762753E-02 ) * va*va*va*va +
        ( 9.66891875E-06 ) * Ta*va*va*va*va +
        ( 2.52785852E-09 ) * Ta*Ta*va*va*va*va +
        ( 4.56306672E-04 ) * va*va*va*va*va +
        ( -1.74202546E-07 ) * Ta*va*va*va*va*va +
        ( -5.91491269E-06 ) * va*va*va*va*va*va +
        ( 3.98374029E-01 ) * DTmrt +
        ( 1.83945314E-04 ) * Ta*DTmrt +
        ( -1.73754510E-04 ) * Ta*Ta*DTmrt +
        ( -7.60781159E-07 ) * Ta*Ta*Ta*DTmrt +
        ( 3.77830287E-08 ) * Ta*Ta*Ta*Ta*DTmrt +
        ( 5.43079673E-10 ) * Ta*Ta*Ta*Ta*Ta*DTmrt +
        ( -2.00518269E-02 ) * va*DTmrt +
        ( 8.92859837E-04 ) * Ta*va*DTmrt +
        ( 3.45433048E-06 ) * Ta*Ta*va*DTmrt +
        ( -3.77925774E-07 ) * Ta*Ta*Ta*va*DTmrt +
        ( -1.69699377E-09 ) * Ta*Ta*Ta*Ta*va*DTmrt +
        ( 1.69992415E-04 ) * va*va*DTmrt +
        ( -4.99204314E-05 ) * Ta*va*va*DTmrt +
        ( 2.47417178E-07 ) * Ta*Ta*va*va*DTmrt +
        ( 1.07596466E-08 ) * Ta*Ta*Ta*va*va*DTmrt +
        ( 8.49242932E-05 ) * va*va*va*DTmrt +
        ( 1.35191328E-06 ) * Ta*va*va*va*DTmrt +
        ( -6.21531254E-09 ) * Ta*Ta*va*va*va*DTmrt +
        ( -4.99410301E-06 ) * va*va*va*va*DTmrt +
        ( -1.89489258E-08 ) * Ta*va*va*va*va*DTmrt +
        ( 8.15300114E-08 ) * va*va*va*va*va*DTmrt +
        ( 7.55043090E-04 ) * DTmrt*DTmrt +
        ( -5.65095215E-05 ) * Ta*DTmrt*DTmrt +
        ( -4.52166564E-07 ) * Ta*Ta*DTmrt*DTmrt +
        ( 2.46688878E-08 ) * Ta*Ta*Ta*DTmrt*DTmrt +
        ( 2.42674348E-10 ) * Ta*Ta*Ta*Ta*DTmrt*DTmrt +
        ( 1.54547250E-04 ) * va*DTmrt*DTmrt +
        ( 5.24110970E-06 ) * Ta*va*DTmrt*DTmrt +
        ( -8.75874982E-08 ) * Ta*Ta*va*DTmrt*DTmrt +
        ( -1.50743064E-09 ) * Ta*Ta*Ta*va*DTmrt*DTmrt +
        ( -1.56236307E-05 ) * va*va*DTmrt*DTmrt +
        ( -1.33895614E-07 ) * Ta*va*va*DTmrt*DTmrt +
        ( 2.49709824E-09 ) * Ta*Ta*va*va*DTmrt*DTmrt +
        ( 6.51711721E-07 ) * va*va*va*DTmrt*DTmrt +
        ( 1.94960053E-09 ) * Ta*va*va*va*DTmrt*DTmrt +
        ( -1.00361113E-08 ) * va*va*va*va*DTmrt*DTmrt +
        ( -1.21206673E-05 ) * DTmrt*DTmrt*DTmrt +
        ( -2.18203660E-07 ) * Ta*DTmrt*DTmrt*DTmrt +
        ( 7.51269482E-09 ) * Ta*Ta*DTmrt*DTmrt*DTmrt +
        ( 9.79063848E-11 ) * Ta*Ta*Ta*DTmrt*DTmrt*DTmrt +
        ( 1.25006734E-06 ) * va*DTmrt*DTmrt*DTmrt +
        ( -1.81584736E-09 ) * Ta*va*DTmrt*DTmrt*DTmrt +
        ( -3.52197671E-10 ) * Ta*Ta*va*DTmrt*DTmrt*DTmrt +
        ( -3.36514630E-08 ) * va*va*DTmrt*DTmrt*DTmrt +
        ( 1.35908359E-10 ) * Ta*va*va*DTmrt*DTmrt*DTmrt +
        ( 4.17032620E-10 ) * va*va*va*DTmrt*DTmrt*DTmrt +
        ( -1.30369025E-09 ) * DTmrt*DTmrt*DTmrt*DTmrt +
        ( 4.13908461E-10 ) * Ta*DTmrt*DTmrt*DTmrt*DTmrt +
        ( 9.22652254E-12 ) * Ta*Ta*DTmrt*DTmrt*DTmrt*DTmrt +
        ( -5.08220384E-09 ) * va*DTmrt*DTmrt*DTmrt*DTmrt +
        ( -2.24730961E-11 ) * Ta*va*DTmrt*DTmrt*DTmrt*DTmrt +
        ( 1.17139133E-10 ) * va*va*DTmrt*DTmrt*DTmrt*DTmrt +
        ( 6.62154879E-10 ) * DTmrt*DTmrt*DTmrt*DTmrt*DTmrt +
        ( 4.03863260E-13 ) * Ta*DTmrt*DTmrt*DTmrt*DTmrt*DTmrt +
        ( 1.95087203E-12 ) * va*DTmrt*DTmrt*DTmrt*DTmrt*DTmrt +
        ( -4.73602469E-12 ) * DTmrt*DTmrt*DTmrt*DTmrt*DTmrt*DTmrt +
        ( 5.12733497E+00 ) * kPa +
        ( -3.12788561E-01 ) * Ta*kPa +
        ( -1.96701861E-02 ) * Ta*Ta*kPa +
        ( 9.99690870E-04 ) * Ta*Ta*Ta*kPa +
        ( 9.51738512E-06 ) * Ta*Ta*Ta*Ta*kPa +
        ( -4.66426341E-07 ) * Ta*Ta*Ta*Ta*Ta*kPa +
        ( 5.48050612E-01 ) * va*kPa +
        ( -3.30552823E-03 ) * Ta*va*kPa +
        ( -1.64119440E-03 ) * Ta*Ta*va*kPa +
        ( -5.16670694E-06 ) * Ta*Ta*Ta*va*kPa +
        ( 9.52692432E-07 ) * Ta*Ta*Ta*Ta*va*kPa +
        ( -4.29223622E-02 ) * va*va*kPa +
        ( 5.00845667E-03 ) * Ta*va*va*kPa +
        ( 1.00601257E-06 ) * Ta*Ta*va*va*kPa +
        ( -1.81748644E-06 ) * Ta*Ta*Ta*va*va*kPa +
        ( -1.25813502E-03 ) * va*va*va*kPa +
        ( -1.79330391E-04 ) * Ta*va*va*va*kPa +
        ( 2.34994441E-06 ) * Ta*Ta*va*va*va*kPa +
        ( 1.29735808E-04 ) * va*va*va*va*kPa +
        ( 1.29064870E-06 ) * Ta*va*va*va*va*kPa +
        ( -2.28558686E-06 ) * va*va*va*va*va*kPa +
        ( -3.69476348E-02 ) * DTmrt*kPa +
        ( 1.62325322E-03 ) * Ta*DTmrt*kPa +
        ( -3.14279680E-05 ) * Ta*Ta*DTmrt*kPa +
        ( 2.59835559E-06 ) * Ta*Ta*Ta*DTmrt*kPa +
        ( -4.77136523E-08 ) * Ta*Ta*Ta*Ta*DTmrt*kPa +
        ( 8.64203390E-03 ) * va*DTmrt*kPa +
        ( -6.87405181E-04 ) * Ta*va*DTmrt*kPa +
        ( -9.13863872E-06 ) * Ta*Ta*va*DTmrt*kPa +
        ( 5.15916806E-07 ) * Ta*Ta*Ta*va*DTmrt*kPa +
        ( -3.59217476E-05 ) * va*va*DTmrt*kPa +
        ( 3.28696511E-05 ) * Ta*va*va*DTmrt*kPa +
        ( -7.10542454E-07 ) * Ta*Ta*va*va*DTmrt*kPa +
        ( -1.24382300E-05 ) * va*va*va*DTmrt*kPa +
        ( -7.38584400E-09 ) * Ta*va*va*va*DTmrt*kPa +
        ( 2.20609296E-07 ) * va*va*va*va*DTmrt*kPa +
        ( -7.32469180E-04 ) * DTmrt*DTmrt*kPa +
        ( -1.87381964E-05 ) * Ta*DTmrt*DTmrt*kPa +
        ( 4.80925239E-06 ) * Ta*Ta*DTmrt*DTmrt*kPa +
        ( -8.75492040E-08 ) * Ta*Ta*Ta*DTmrt*DTmrt*kPa +
        ( 2.77862930E-05 ) * va*DTmrt*DTmrt*kPa +
        ( -5.06004592E-06 ) * Ta*va*DTmrt*DTmrt*kPa +
        ( 1.14325367E-07 ) * Ta*Ta*va*DTmrt*DTmrt*kPa +
        ( 2.53016723E-06 ) * va*va*DTmrt*DTmrt*kPa +
        ( -1.72857035E-08 ) * Ta*va*va*DTmrt*DTmrt*kPa +
        ( -3.95079398E-08 ) * va*va*va*DTmrt*DTmrt*kPa +
        ( -3.59413173E-07 ) * DTmrt*DTmrt*DTmrt*kPa +
        ( 7.04388046E-07 ) * Ta*DTmrt*DTmrt*DTmrt*kPa +
        ( -1.89309167E-08 ) * Ta*Ta*DTmrt*DTmrt*DTmrt*kPa +
        ( -4.79768731E-07 ) * va*DTmrt*DTmrt*DTmrt*kPa +
        ( 7.96079978E-09 ) * Ta*va*DTmrt*DTmrt*DTmrt*kPa +
        ( 1.62897058E-09 ) * va*va*DTmrt*DTmrt*DTmrt*kPa +
        ( 3.94367674E-08 ) * DTmrt*DTmrt*DTmrt*DTmrt*kPa +
        ( -1.18566247E-09 ) * Ta*DTmrt*DTmrt*DTmrt*DTmrt*kPa +
        ( 3.34678041E-10 ) * va*DTmrt*DTmrt*DTmrt*DTmrt*kPa +
        ( -1.15606447E-10 ) * DTmrt*DTmrt*DTmrt*DTmrt*DTmrt*kPa +
        ( -2.80626406E+00 ) * kPa*kPa +
        ( 5.48712484E-01 ) * Ta*kPa*kPa +
        ( -3.99428410E-03 ) * Ta*Ta*kPa*kPa +
        ( -9.54009191E-04 ) * Ta*Ta*Ta*kPa*kPa +
        ( 1.93090978E-05 ) * Ta*Ta*Ta*Ta*kPa*kPa +
        ( -3.08806365E-01 ) * va*kPa*kPa +
        ( 1.16952364E-02 ) * Ta*va*kPa*kPa +
        ( 4.95271903E-04 ) * Ta*Ta*va*kPa*kPa +
        ( -1.90710882E-05 ) * Ta*Ta*Ta*va*kPa*kPa +
        ( 2.10787756E-03 ) * va*va*kPa*kPa +
        ( -6.98445738E-04 ) * Ta*va*va*kPa*kPa +
        ( 2.30109073E-05 ) * Ta*Ta*va*va*kPa*kPa +
        ( 4.17856590E-04 ) * va*va*va*kPa*kPa +
        ( -1.27043871E-05 ) * Ta*va*va*va*kPa*kPa +
        ( -3.04620472E-06 ) * va*va*va*va*kPa*kPa +
        ( 5.14507424E-02 ) * DTmrt*kPa*kPa +
        ( -4.32510997E-03 ) * Ta*DTmrt*kPa*kPa +
        ( 8.99281156E-05 ) * Ta*Ta*DTmrt*kPa*kPa +
        ( -7.14663943E-07 ) * Ta*Ta*Ta*DTmrt*kPa*kPa +
        ( -2.66016305E-04 ) * va*DTmrt*kPa*kPa +
        ( 2.63789586E-04 ) * Ta*va*DTmrt*kPa*kPa +
        ( -7.01199003E-06 ) * Ta*Ta*va*DTmrt*kPa*kPa +
        ( -1.06823306E-04 ) * va*va*DTmrt*kPa*kPa +
        ( 3.61341136E-06 ) * Ta*va*va*DTmrt*kPa*kPa +
        ( 2.29748967E-07 ) * va*va*va*DTmrt*kPa*kPa +
        ( 3.04788893E-04 ) * DTmrt*DTmrt*kPa*kPa +
        ( -6.42070836E-05 ) * Ta*DTmrt*DTmrt*kPa*kPa +
        ( 1.16257971E-06 ) * Ta*Ta*DTmrt*DTmrt*kPa*kPa +
        ( 7.68023384E-06 ) * va*DTmrt*DTmrt*kPa*kPa +
        ( -5.47446896E-07 ) * Ta*va*DTmrt*DTmrt*kPa*kPa +
        ( -3.59937910E-08 ) * va*va*DTmrt*DTmrt*kPa*kPa +
        ( -4.36497725E-06 ) * DTmrt*DTmrt*DTmrt*kPa*kPa +
        ( 1.68737969E-07 ) * Ta*DTmrt*DTmrt*DTmrt*kPa*kPa +
        ( 2.67489271E-08 ) * va*DTmrt*DTmrt*DTmrt*kPa*kPa +
        ( 3.23926897E-09 ) * DTmrt*DTmrt*DTmrt*DTmrt*kPa*kPa +
        ( -3.53874123E-02 ) * kPa*kPa*kPa +
        ( -2.21201190E-01 ) * Ta*kPa*kPa*kPa +
        ( 1.55126038E-02 ) * Ta*Ta*kPa*kPa*kPa +
        ( -2.63917279E-04 ) * Ta*Ta*Ta*kPa*kPa*kPa +
        ( 4.53433455E-02 ) * va*kPa*kPa*kPa +
        ( -4.32943862E-03 ) * Ta*va*kPa*kPa*kPa +
        ( 1.45389826E-04 ) * Ta*Ta*va*kPa*kPa*kPa +
        ( 2.17508610E-04 ) * va*va*kPa*kPa*kPa +
        ( -6.66724702E-05 ) * Ta*va*va*kPa*kPa*kPa +
        ( 3.33217140E-05 ) * va*va*va*kPa*kPa*kPa +
        ( -2.26921615E-03 ) * DTmrt*kPa*kPa*kPa +
        ( 3.80261982E-04 ) * Ta*DTmrt*kPa*kPa*kPa +
        ( -5.45314314E-09 ) * Ta*Ta*DTmrt*kPa*kPa*kPa +
        ( -7.96355448E-04 ) * va*DTmrt*kPa*kPa*kPa +
        ( 2.53458034E-05 ) * Ta*va*DTmrt*kPa*kPa*kPa +
        ( -6.31223658E-06 ) * va*va*DTmrt*kPa*kPa*kPa +
        ( 3.02122035E-04 ) * DTmrt*DTmrt*kPa*kPa*kPa +
        ( -4.77403547E-06 ) * Ta*DTmrt*DTmrt*kPa*kPa*kPa +
        ( 1.73825715E-06 ) * va*DTmrt*DTmrt*kPa*kPa*kPa +
        ( -4.09087898E-07 ) * DTmrt*DTmrt*DTmrt*kPa*kPa*kPa +
        ( 6.14155345E-01 ) * kPa*kPa*kPa*kPa +
        ( -6.16755931E-02 ) * Ta*kPa*kPa*kPa*kPa +
        ( 1.33374846E-03 ) * Ta*Ta*kPa*kPa*kPa*kPa +
        ( 3.55375387E-03 ) * va*kPa*kPa*kPa*kPa +
        ( -5.13027851E-04 ) * Ta*va*kPa*kPa*kPa*kPa +
        ( 1.02449757E-04 ) * va*va*kPa*kPa*kPa*kPa +
        ( -1.48526421E-03 ) * DTmrt*kPa*kPa*kPa*kPa +
        ( -4.11469183E-05 ) * Ta*DTmrt*kPa*kPa*kPa*kPa +
        ( -6.80434415E-06 ) * va*DTmrt*kPa*kPa*kPa*kPa +
        ( -9.77675906E-06 ) * DTmrt*DTmrt*kPa*kPa*kPa*kPa +
        ( 8.82773108E-02 ) * kPa*kPa*kPa*kPa*kPa +
        ( -3.01859306E-03 ) * Ta*kPa*kPa*kPa*kPa*kPa +
        ( 1.04452989E-03 ) * va*kPa*kPa*kPa*kPa*kPa +
        ( 2.47090539E-04 ) * DTmrt*kPa*kPa*kPa*kPa*kPa +
        ( 1.48348065E-03 ) * kPa*kPa*kPa*kPa*kPa*kPa;
}



Foam::scalar Foam::comfortFunctions::correctSET
(
    const scalar& pressure,
    const scalar& relativeHumidity,

    scalar& velocity,
    scalar& temperature,
    scalar& radiantTemperature,
    const bool occupantControl
 )
{
    scalar ASHSET = -1;

    if (velocity > 0.2)
    {
            //Adjust velocity if occupant controlled
        if (!occupantControl)
        {
            scalar to = TOp
                (
                    velocity,
                    temperature,
                    radiantTemperature
                 );

            if (to < 23.0)
            {
                velocity = 0.2;
            }
            else if (to <= 25.5)
            {
                velocity = 50.49-(4.4047*to)+0.096425*sqr(to);
            }
            else
            {
                velocity = 0.8;
            }
        }

        ASHSET = SET
            (
                velocity,
                temperature,
                pressure,
                radiantTemperature,
                relativeHumidity
             );

            //reset velocity
        velocity = 0.1;
        scalar Tx = temperature;
        scalar rTx = radiantTemperature;

        scalar SETx = GREAT;
        scalar tDrop = 0.0;
        scalar tDropFactor = 5.0;

        label maxIter = 100;
        label iter = 0;

        while (SETx > ASHSET && iter++ < maxIter)
        {
            tDrop += tDropFactor;
            Tx = temperature - tDrop;
            rTx = radiantTemperature - tDrop;
            SETx = SET
                (
                    velocity,
                    Tx,
                    pressure,
                    rTx,
                    relativeHumidity
                 );
        }

            //perform bi-section to get temperature drop
        scalar tol = 0.05;

        scalar dTUpper = tDrop-tDropFactor;
        scalar dTLower = tDrop;

        iter = 0;
        while (mag(SETx-ASHSET) > tol && iter++ < maxIter)
        {
            tDrop = 0.5*(dTUpper + dTLower);
            Tx = temperature - tDrop;
            rTx = radiantTemperature - tDrop;
            SETx = SET
                (
                    velocity,
                    Tx,
                    pressure,
                    rTx,
                    relativeHumidity
                 );

            if (SETx < ASHSET)
            {
                dTLower = tDrop;
            }
            else
            {
                dTUpper = tDrop;
            }
        }

            //Update temperature and relative temperature
            // for PMV, PPD corrected calculation
        temperature = Tx;
        radiantTemperature = rTx;
    }
    else
    {
            // use uncorrected equivalent temperature
        ASHSET = SET
            (
                velocity,
                temperature,
                pressure,
                radiantTemperature,
                relativeHumidity
             );
    }

    return ASHSET;
}


Foam::scalar Foam::comfortFunctions::SET
(
    const scalar& velocity,
    const scalar& temperature,
    const scalar& pressure,
    const scalar& radiantTemperature,
    const scalar& relativeHumidity
 )
{
        //not valid above 50 C
    if (temperature > 323.15)
    {
        return 50;
    }

    scalar clo1 = 0.25;
    scalar weight = 69.9; //kg
    scalar surfaceArea = 1.8258; //m2
    scalar metConvert = 58.15; //W/m2

    scalar csw = 170.0;
    scalar cdil = 120.0;
    scalar cstr = 0.5;
    label ltime = 60;

    scalar met = met_/metConvert;
    scalar clo = clo_/0.155;

    scalar skinNeutralTemp = 33.7;
    scalar coreNeutralTemp = 36.8;
    scalar bodyNeutralTemp = 36.49;
    scalar skinBloodNeutral = 6.3;

    scalar rcl = 0.155 * clo;
    scalar facl = 1.0 + 0.15 * clo;

    scalar icl(clo <= 0 ? 1.0 : 0.45);

        //Convert temperature to celsius
    scalar TA = temperature - 273.15;
        //Convert radiant temperature to celsius
    scalar TR =  radiantTemperature - 273.15;
    scalar airVelocity = max(velocity, 0.1);
    scalar vaporP =  relativeHumidity * svpTorr(temperature)/100.0;

    scalar atmP = pressure * 9.86923e-6;

    scalar skinTemp = skinNeutralTemp; //Initial values
    scalar coreTemp = coreNeutralTemp;
    scalar skinBloodFlow = skinBloodNeutral;

    scalar mshiv = 0.0;
    scalar alfa = 0.1;

    scalar rm = met * metConvert;
    scalar mw = met * metConvert;

    scalar lr = 2.2/atmP;
    scalar wcrit
        (
            clo <= 0
            ? 0.38 * pow(airVelocity, -0.29)
            : 0.59 * pow(airVelocity, -0.08)
         );

    scalar esk = 0.1 * met;
    scalar chc = 3.0 * pow(atmP, 0.53);
    scalar chcV = 8.600001 * pow((airVelocity * atmP), 0.53);
    chc = max(chc, chcV);
    scalar chr = 4.7;
    scalar ctc = chr + chc;
    scalar ra = 1.0/(facl * ctc);
    scalar top = (chr * TR + chc * TA)/ctc;
    scalar tcl = top + (skinTemp - top)/(ctc * (ra + rcl));

    scalar tcl_old = tcl;
    bool flag = true;
    scalar dry,hfcs,eres,cres,scr,ssk,tcsk,tccr,dtsk,dtcr,tb;
    scalar sksig,warms,colds,crsig,warmc,coldc,bdsig,warmb;
    scalar regsw,ersw,rea,recl,emax,prsw,pwet,edif;

    for (int tim = 1; tim <= ltime; tim++)
    {
        do
        {
            if (flag)
            {
                tcl_old = tcl;
                chr = 4.0 * sigma_
                    * pow(((tcl + TR)/2.0 + 273.15), 3.0) * 0.72;
                ctc = chr + chc;
                ra = 1.0/(facl * ctc);
                top = (chr * TR + chc * TA)/ctc;
            }
            tcl = (ra * skinTemp + rcl * top)/(ra + rcl);
            flag = true;
        } while (mag(tcl - tcl_old) > 0.01);
        flag = false;

        dry = (skinTemp - top)/(ra + rcl);
        hfcs = (coreTemp - skinTemp) * (5.28 + 1.163 * skinBloodFlow);
        eres = 0.0023 * mw * (44.0 - vaporP);
        cres = 0.0014 * mw * (34.0 - TA);
        scr = mw - hfcs - eres - cres - work_;
        ssk = hfcs - dry - esk;
        tcsk = 0.97 * alfa * weight;
        tccr = 0.97 * (1 - alfa) * weight;
        dtsk = (ssk * surfaceArea)/(tcsk * 60.0);
        dtcr = scr * surfaceArea/(tccr * 60.0);
        skinTemp = skinTemp + dtsk;
        coreTemp = coreTemp + dtcr;
        tb = alfa * skinTemp + (1 - alfa) * coreTemp;
        sksig = skinTemp - skinNeutralTemp;

        if (sksig > 0)
        {
            warms = sksig;
            colds = scalar(0.0);
        }
        else
        {
            warms = 0.0;
            colds = scalar(-1.0) * sksig;
        }
        crsig = (coreTemp - coreNeutralTemp);
        if (crsig > 0)
        {
            warmc = crsig;
            coldc = scalar(0.0);
        }
        else
        {
            warmc = scalar(0.0);
            coldc = scalar(-1.0) * crsig;
        }
        bdsig = tb - bodyNeutralTemp;

        warmb = (bdsig > 0) * bdsig;
        skinBloodFlow = (skinBloodNeutral + cdil *warmc)
            /(1 + cstr * colds);
        skinBloodFlow =
            max(scalar(0.5), min(scalar(90.0), skinBloodFlow));
        regsw = csw *warmb * exp(warms/10.7);
        regsw = min(regsw, scalar(500.0));
        ersw = 0.68 * regsw;
        rea = 1.0/(lr * facl * chc);
        recl = rcl/(lr * icl);
        emax = (svpTorr(skinTemp+273.15) - vaporP)
            /(rea + recl);
        prsw = ersw/emax;
        pwet = 0.06 + 0.94 * prsw;
        edif = pwet * emax - ersw;
        esk = ersw + edif;
        if (pwet > wcrit)
        {
            pwet = wcrit;
            prsw = wcrit/0.94;
            ersw = prsw * emax;
            edif = 0.06 * (1.0 - prsw) * emax;
            esk = ersw + edif;
        }

        if (emax < 0)
        {
            edif = 0;
            ersw = 0;
            pwet = wcrit;
            prsw = wcrit;
            esk = emax;
        }
        esk = ersw + edif;
        mshiv = 19.4 * colds * coldc;
        mw = rm + mshiv;
        alfa = 0.0417737 + 0.7451833/(skinBloodFlow + 0.585417);
    }

    scalar hsk = dry + esk;
    scalar rn = mw - work_;
    scalar ecomf = 0.42 * (rn - (1 * metConvert));
    if (ecomf < 0.0) ecomf = 0.0;
    emax = emax * wcrit;
    scalar w = pwet;
    scalar pssk = svpTorr(skinTemp+273.15);
    scalar chrS = chr;

    scalar chcS;

    if (met < 0.85)
    {
        chcS = 3.0;
    }
    else
    {
        chcS = 5.66 * pow((met - 0.85), 0.39);
        chcS = max(chcS, scalar(3.0));
    }
    scalar ctcS = chcS + chrS;
    scalar rclos = 1.52
        /((met - work_/metConvert) + 0.6944) - 0.1835;
    scalar rcls = 0.155 * rclos;
    scalar facls = 1.0 + clo1 * rclos;
    scalar fcls = 1.0/(1.0 + 0.155 * facls * ctcS * rclos);
    scalar ims = 0.45;
    scalar icls = ims * chcS/ctcS
        * (1 - fcls)/(chcS/ctcS - fcls * ims);
    scalar raS = 1.0/(facls * ctcS);
    scalar reaS = 1.0/(lr * facls * chcS);
    scalar reclS = rcls/(lr * icls);
    scalar hd_s = 1.0/(raS + rcls);
    scalar he_s = 1.0/(reaS + reclS);

    scalar tol = .0001;
    scalar dx = 100.0;
    scalar seti, err1, err2;
    scalar set_old = skinTemp - hsk/hd_s;
    while (mag(dx) > .01)
    {
        err1 =
            (
                hsk - hd_s * (skinTemp - set_old) -w * he_s
                * (pssk - 0.5 *svpTorr(set_old+273.15))
             );
        err2 =
            (
                hsk - hd_s * (skinTemp - (set_old + tol)) -w * he_s
                * (pssk - 0.5 *svpTorr((set_old+273.15 + tol)))
             );
        seti = set_old - tol * err1/(err2 - err1);
        dx = seti - set_old;

        set_old = seti;
    }

        //limit maximum equivalent temperature
    return min(seti, 50);
}


Foam::scalarField Foam::comfortFunctions::PHSiso7933
(
    const scalar& va,
    const scalar& temperature,
    const scalar& radiantTemperature,
    const scalar& relativeHumidity,
    const scalar& weight,
    const scalar& height,
    const Switch& drink,
    const scalar& duration,
    const scalar& posture,
    const Switch& accl,
    const scalar& thetaW,
    const scalar& walkSpeed
 )
{
    scalar Ta = temperature - 273.15;
    scalar Tmrt = radiantTemperature - 273.15;

        //- Saturation Vapor Pressure - Arden Buck equation
    scalar Pvs = 0.61121*Foam::exp((18.678 - Ta/234.5)*(Ta/(Ta + 257.14)));
    scalar Pa = Pvs*relativeHumidity/100;
        //scalar Pa = relativeHumidity/1000/(0.622+relativeHumidity/1000)*101325/1000;

    scalar clo = clo_ / 0.155;

    scalar Adu = 0.202*(Foam::pow(weight,0.425))*(Foam::pow(height,0.725));
    scalar Ardu = 0.0;
    scalar spHeat = 57.83*weight/Adu;
    scalar SWp = 0.0;
    scalar Tre = 36.8;
    scalar Tcr = 36.8;
    scalar TcrEq = Tcr;
    scalar TcrEqm = Tcr;
    scalar Tsk = 34.1;
    scalar TskTcrWg = 0.3;
    scalar Dmax50 = 0.075*weight*1000;
    scalar Dmax95 = 0.05*weight*1000;
    scalar ConstTeq = Foam::exp(-1.0/10.0);
    scalar ConstTsk = Foam::exp(-1.0/3.0);
    scalar ConstSW = Foam::exp(-1.0/10.0);
    scalar imst = 0.38;
    scalar Ap = 0.54;
    scalar Fr = 0.97;
    Switch defDir;
    Switch defSpeed;

    if (thetaW != 0.0)
    {
        defDir = true;
    } else
    {
        defDir = false;
    }

    if (walkSpeed == 0.0)
    {
        defSpeed = false;
    } else
    {
        defSpeed = true;
    }

    scalar Texp,Cres,Eres,Ep,k,wp,Rtdyn,Emax,Ereq,Itotdyn,imdyn;
    scalar Var = 0.0;
    scalar walkAux = 0.0;
    scalar wSpeed = walkSpeed;
    scalar dStoreq,Conv,Rad,hcDyn,fcl,Tcl,auxR,FclR,Icldyn;
    scalar TskEq,TskEqCl,TskEqNu,Psk;

    scalar timePHS(0);
    scalar Tsk_0(0);
    scalar Tre_0(0);
    scalar Tcr_0(0);
    scalar TcrEq_0(0);
    scalar TskTcrWg_0(0);
    scalar SWtot(0);
    scalar Hr(0);

    scalar DlimTre(0);
    scalar DlimLoss50(0);
    scalar DlimLoss95(0);
    scalar SWtotg(0);

    for (timePHS = 1; timePHS < duration+1; timePHS++)
    {
        Tsk_0 = Tsk;
        Tre_0 = Tre;
        Tcr_0 = Tcr;
        TcrEq_0 = TcrEq;
        TskTcrWg_0 = TskTcrWg;

            //- Effective radiating area of the body
        if (posture == 1)
        {
            Ardu = 0.7;
        }
        if (posture == 2)
        {
            Ardu = 0.77;
        }
        if (posture == 3)
        {
            Ardu = 0.67;
        }

            //- Evaluation of the maximum sweet rate as a function of met
        scalar SWmax = (met_ - 32.0)*Adu;
        if (SWmax > 400.0)
        {
            SWmax = 400.0;
        }

        if (SWmax < 250.0)
        {
            SWmax = 250.0;
        }

            //- For acclimatised subject, the maximum sweet rate is greater by 25%
        scalar Wmax = 0.0;
        if (accl == true)
        {
            SWmax = SWmax*1.25;
            Wmax = 1.0;
        }
        else
        {
            Wmax = 0.85;
        }

            //- Equilibrium core temperature as function of met
        TcrEqm = 0.0036*(met_-55.0) + 36.8;

            //- Core temperature at this minute, by exponential averaging
        TcrEq = TcrEq_0*ConstTeq + TcrEqm*(1-ConstTeq);

            //- Heat storage associated with this core temperature increase
            // during the last minute
        dStoreq = spHeat*(TcrEq - TcrEq_0)*(1-TskTcrWg_0);

            //- Skin temperature in equilibrium
        TskEqCl = 12.165 + 0.02017*Ta + 0.04361*Tmrt
            + 0.19354*Pa - 0.25315*va;
        TskEqCl = TskEqCl + 0.005346*met_ + 0.51274*Tre; // clothed model
        TskEqNu = 7.191 + 0.064*Ta + 0.061*Tmrt + 0.198*Pa - 0.348*va;
        TskEqNu = TskEqNu + 0.616*Tre; // nude model

            //- Skin temperature at this minute, as a function of
            // clothing insulation
        if (clo >= 0.6)
        {
            TskEq = TskEqCl;
        }
        else if (clo <= 0.2)
        {
            TskEq = TskEqNu;
        }
        else
        {
            TskEq = TskEqNu + 2.5*(TskEqCl - TskEqNu)*(clo - 0.2);
        }

            //- Skin temperature at this minute, by exponential averaging
        Tsk = Tsk_0*ConstTsk + TskEq*(1-ConstTsk);
            //- Saturated water vapour pressure at the surface of the skin
        Psk = 0.6105*Foam::exp(17.27*Tsk/(Tsk+237.3));

            //- Clothing influence on exchange coefficients
        scalar Iclst = clo*0.155; //static clothing insulation
        fcl = 1 + 0.3*clo; //clothing area factor
        scalar Iast = 0.111;       //thermal insulation in quiet air
        scalar Itotst = Iclst + Iast/fcl;

            //- Relative velocities due to air velocity and movements
        if (defSpeed == true)
        {
            if (defDir == true)
            {
                Var = fabs(va - wSpeed*cos(M_PI*thetaW/180.0));
            } else
            {
                if (va < wSpeed)
                {
                    Var = wSpeed;
                } else
                {
                    Var = va;
                }
            }
        }
        else
        {
            wSpeed = 0.0052*(met_ - 58.0);
            if (wSpeed > 0.7)
            {
                wSpeed = 0.7;
            }
            Var = va;
        }

            //- Dynamic Clothing Insulation: correction for wind (Var)
            // and walking (walkSpeed)
        scalar vaux = Var;
        walkAux = wSpeed;
        if (Var > 3.0)
        {
            vaux = 3.0;
        }
        if (wSpeed > 1.5)
        {
            walkAux = 1.5;
        }
        scalar CORcl = 1.044*Foam::exp((0.066*vaux - 0.398)*vaux
                                       + (0.094*walkAux - 0.378)*walkAux);
        if (CORcl > 1.0)
        {
            CORcl = 1.0;
        }
        scalar CORia = Foam::exp((0.047*Var - 0.472)*Var
                                 + (0.117*walkAux - 0.342)*walkAux);
        if (CORia > 1.0)
        {
            CORia = 1.0;
        }
        scalar CORtot = CORcl;
        if (clo <= 0.6)
        {
            CORtot = ((0.6 - clo)*CORia + clo*CORcl)/0.6;
        }

        Itotdyn = Itotst*CORtot;
        scalar IAdyn = CORia*Iast;

        Icldyn = Itotdyn - IAdyn/fcl;

            //- Permeability index: correction for wind and walking
        scalar CORe = (2.6*CORtot - 6.5)*CORtot + 4.9;
        imdyn = imst*CORe;
        if (imdyn > 0.9)
        {
            imdyn = 0.9;
        }

        Rtdyn = Itotdyn/imdyn/16.7;

            //- Heat Exchanges
        Texp = 28.56 + 0.115*Ta + 0.641*Pa;
        Cres = 0.001516 * met_ * (Texp - Ta);
        Eres = 0.00127 * met_ * (59.34 + 0.53*Ta - 11.63*Pa);

            //- Calculation of Mean Temperature of the Clothing
        scalar Z = 3.5 + 5.2*Var;
        if (Var > 1.0)
        {
            Z = 8.7*Foam::pow(Var,0.6);
        }
        hcDyn = 2.38*Foam::pow(fabs(Tsk - Ta),0.25);
        if (Z > hcDyn)
        {
            hcDyn = Z;
        }

        auxR = 0.0000000567*Ardu;
        FclR = (1 - Ap)*0.97 + Ap*Fr;
        Tcl = Tmrt + 0.1;

            //- Tcl iteration
        scalar Tcl1 = 0;

        while (fabs(Tcl-Tcl1) > 0.001)
        {
            Hr = FclR*auxR*(Foam::pow((Tcl + 273.0),4)
                            - Foam::pow((Tmrt + 273.0),4))/(Tcl - Tmrt);
            Tcl1 = ((fcl*(hcDyn*Ta + Hr*Tmrt) + Tsk/Icldyn))
                /(fcl*(hcDyn + Hr) + 1.0/Icldyn);
            Tcl = (Tcl + Tcl1)/2;
        }

        Conv = fcl*hcDyn*(Tcl - Ta);
        Rad = fcl*Hr*(Tcl - Tmrt);
            //- maximum evaporate rate
        Emax = (Psk - Pa)/Rtdyn;
            //- required evaporate rate
        Ereq = met_ - dStoreq - work_ - Cres - Eres - Conv - Rad;

            //- required wettedness
        scalar wreq = Ereq/Emax;
        scalar SWreq = 0.0;
            //-required sweat rate
        if (Ereq <= 0.0)
        {
            Ereq = 0.0;
            SWreq = 0.0;
        }
        else if (Emax <= 0.0)
        {
            Emax = 0.0;
            SWreq = SWmax;
        }
        else if (wreq >= 1.7)
        {
            wreq = 1.7;
            SWreq = SWmax;
        }
        else
        {
            scalar Eveff = (1 - Foam::pow(wreq,2)/2); //here
            if (wreq > 1)
            {
                Eveff = Foam::pow((2 - wreq),2)/2;
            }
            SWreq = Ereq / Eveff; //here
            if (SWreq > SWmax)
            {
                SWreq = SWmax;
            }
        }

            //- Predicted Sweat Rate by exponential averaging
        SWp = SWp*ConstSW + SWreq*(1-ConstSW);

        if (SWp <= 0.0)
        {
            Ep = 0.0;
            SWp = 0.0;
        }
        else
        {
            k = Emax/SWp;
            wp = 1;
            if (k >= 0.5)
            {
                wp = -k + sqrt(k*k+2);
            }
            if (wp > Wmax)
            {
                wp = Wmax;
            }
            Ep = wp*Emax;
        }

            //- Heat Storage
        scalar dStorage = Ereq - Ep + dStoreq;

            //- Prediction of the core temperature
        scalar Tcr1 = Tcr_0;
        scalar Tcr_temp;

        do
        {
            TskTcrWg = 0.3 - 0.09*(Tcr1 - 36.8);
            if (TskTcrWg > 0.3)
            {
                TskTcrWg = 0.3;
            }
            if (TskTcrWg < 0.1)
            {
                TskTcrWg = 0.1;
            }

            Tcr = dStorage/spHeat + Tsk_0*TskTcrWg_0/2.0 - Tsk*TskTcrWg/2.0;
            Tcr = (Tcr + Tcr_0*(1.0 - TskTcrWg_0/2.0))/(1.0-TskTcrWg/2.0);
            Tcr_temp = Tcr1;
            Tcr1 = (Tcr1 + Tcr)/2.0;
        }while(mag(Tcr - Tcr_temp) > 0.001);

            //- Prediction of the rectal temperature
        Tre = Tre_0 + (2.0*Tcr - 1.962*Tre_0 - 1.31)/9.0;

        if (DlimTre == 0 && Tre >= 38.0)
        {
            DlimTre = timePHS;
        }

            //- Total Water Loss rate during the minute
        SWtot = SWtot + SWp + Eres;
        SWtotg = SWtot*2.67*Adu/1.8/60.0;
        if (DlimLoss50 == 0 && SWtotg >= Dmax50)
        {
            DlimLoss50 = timePHS;
        }

        if (DlimLoss95 == 0 && SWtotg >= Dmax95)
        {
            DlimLoss95 = timePHS;
        }

        if (drink == false)
        {
            DlimLoss95 = DlimLoss95 * 0.6;
            DlimLoss50 = DlimLoss95;
        }
    }
        //- end of loop on duration

        //- Dlim computation
    if (DlimLoss50 == 0)
    {
        DlimLoss50 = duration;
    }
    if (DlimLoss95 == 0)
    {
        DlimLoss95 = duration;
    }
    if (DlimTre == 0)
    {
        DlimTre = duration;
    }

    scalarField iso7933Metrics(5);

        //- return iso-7933 cell metrics
    iso7933Metrics[0] = Tre;
    iso7933Metrics[1] = SWtotg;
    iso7933Metrics[2] = DlimTre;
    iso7933Metrics[3] = DlimLoss50;
    iso7933Metrics[4] = DlimLoss95;

    return iso7933Metrics;
}

Foam::scalar Foam::comfortFunctions::fTnwb
(
    scalar& Tit,
    const scalar& pressure,
    const scalar& temperature,
    const scalar& relativeHumidity,
    const scalar& va,
    const scalar& rho,
    const scalar& nu,
    const scalar& Cp,
    const scalar& kappa,
    const scalar& SR,
    const scalar& zenithAngle
)
{
    scalar Tk = temperature;
    scalar RH = relativeHumidity;
    scalar hectoPressure = pressure/100.0;

    // Constants
    scalar MH2O = 18.015;
    scalar Mair = 28.97;

        // Wick constants
    scalar emissWick = 0.95;
    scalar albWick = 0.4;
    scalar lengthWick = 0.0254;
    scalar wickDiam = 0.007;

        // Surface constants
    scalar emissSfc = 0.999;
    scalar albSfc = 0.4;


    // Saturation Vapor Pressure - modified Arden Buck equation
    scalar Pvs = 1.004*6.1121*Foam::exp(17.502*(Tk - 273.15)/(Tk - 32.18));
    scalar partialPa = RH*Pvs;
    scalar partialPw = 6.1121*Foam::exp(17.502*(Tit - 273.15)/(Tit - 32.18));

    // Atmospheric emissivity
    scalar emissAtm = 0.575*pow(RH*Pvs,scalar(0.143));

    // Constants for diffusivity of water vapor in air (BSL, page 505.)

    scalar pcrit13 = pow((36.4*218.0),(1.0/3.0));
    scalar tcrit512 = pow((132.0*647.3),(5.0/12.0));
    scalar Tcrit12 = pow((132.0*647.3),0.5);
    scalar Mmix = pow((1.0/28.97 + 1.0/18.015),0.5);
    scalar Tref = 0.5*(Tit + Tk);

    scalar diffWaterVapor =
    (
        0.000364*pow((Tref/Tcrit12),2.334)
       *pcrit13*tcrit512*Mmix/(hectoPressure/1013.25)*0.0001
    );

    // Radiation Parameters

    scalar tanZA = tan(zenithAngle);
    scalar cosZA = cos(zenithAngle);

        // fraction of the total horizontal solar irradiance SR
    scalar fDir = 0.0;

    if (zenithAngle <= 1.56) // 89.5 degrees (1.56 rad)
    {
        // assumed mean earth-sun distance = 1 A.U
        scalar SRstar = max(SR,0.000001)/(1367.0*cosZA);
        fDir = exp(3 - 1.34*SRstar - 1.65/SRstar);
    }

    // Heat transfer coefficient of cylinder in air

    scalar Pr = nu*rho*Cp/kappa;
    scalar Reynolds = va*wickDiam/nu;
    scalar Nusselt = 0.281*pow(Reynolds,scalar(0.6))*pow(Pr,scalar(0.44));

    scalar htc = Nusselt*kappa/wickDiam;

    // Heat of evaporation
    scalar deltaH = ((313.15 - Tit)/30.0)*(-71100.0) + 2407300.0;

    // Radiative Heating term
    scalar deltaFatm =
    (
        sigma_*emissWick*(0.5*(emissAtm*pow(Tk,4.0)
       +emissSfc*pow(Tk,4.0)) - pow(Tit,4.0))
       +(1-albWick)*SR*((1-fDir)*(1+0.25*wickDiam/lengthWick)
       +((tanZA)/3.1416 + 0.25*wickDiam/lengthWick)*fDir + albSfc)
    );

    scalar Sc = nu/diffWaterVapor;

    return
    (
        Tk - (deltaH/Cp)*(MH2O/Mair)*(Pr/Sc)
       *((partialPw - partialPa)/(hectoPressure - partialPw))
       +deltaFatm/htc
    );
}

Foam::scalar Foam::comfortFunctions::fTg
(
    scalar& Tit,
    const scalar& temperature,
    const scalar& relativeHumidity,
    const scalar& va,
    const scalar& rho,
    const scalar& nu,
    const scalar& Cp,
    const scalar& kappa,
    const scalar& SR,
    const scalar& zenithAngle
)
{
    scalar Tk = temperature;
    scalar RH = relativeHumidity;

    // Constants

        // Globe constants

    scalar emissGlobe = 0.95;
    scalar albGlobe = 0.05;
    scalar globeDiam = 0.05;
    scalar areaGlobe = 3.1416*pow(globeDiam,2);

        // Surface constants

    //scalar emissSfc = 0.999;
    scalar albSfc = 0.4;


    //- Saturation Vapor Pressure - modified Arden Buck equation
    scalar Pvs = 1.004*6.1121*Foam::exp(17.502*(Tk - 273.15)/(Tk - 32.18));

    // Atmospheric emissivity
    scalar emissAtm = 0.575*pow(RH*Pvs,scalar(0.143));


    // Radiation Parameters

    scalar cosZA = cos(zenithAngle);

        // fraction of the total horizontal solar irradiance SR
    scalar fDir = 0.0;

    if (zenithAngle <= 1.56) // 89.5 degrees (1.56 rad)
    {
        // assumed mean earth-sun distance = 1 A.U
        scalar SRstar = max(SR,0.000001)/(1367.0*cosZA);
        fDir = exp(3.0 - 1.34*SRstar - 1.65/SRstar);
    }

    // HTC of sphere on air
    scalar Pr = nu*rho*Cp/kappa;
    scalar Reynolds = va*globeDiam/nu;
    scalar Nusselt = 2 + 0.6*pow(Reynolds,scalar(0.5))*pow(Pr,scalar(0.3333));

    scalar htc = Nusselt*kappa/globeDiam;


    // heat loss by globe surface due to radiation
    scalar qrG = areaGlobe * emissGlobe * sigma_ * pow(Tit,4);

    // heat loss by globe surface due to convection
    scalar qhG = areaGlobe * htc * (Tit - Tk);

    // heat gained by globe surface due to radiation emmited by atmosphere
    scalar qrAtm = scalar(0.5) * areaGlobe * emissAtm * sigma_ * emissGlobe
        * pow(Tk,4);

    // heat gained by globe surface due to radiation emmited by ground surface
    scalar qrSfc = scalar(0.5) * areaGlobe * sigma_ * emissGlobe * pow(Tk,4);

    // heat gained by globe surface due to radiation emmited by sun (diffuse parcel)
    scalar qrSunDiff = scalar(0.5) * areaGlobe * (scalar(1.0) - albGlobe)
        * (1 - fDir) * SR;

    // heat gained by globe surface due to radiation emmited by sun (direct parcel)
    scalar qrSunDir = scalar(0.25) * areaGlobe * (scalar(1.0) - albGlobe)
        * fDir * SR / cosZA;

    // heat gained by globe surface due to radiation emmited by sun (reflected parcel)
    scalar qrSunRef = scalar(0.5) * areaGlobe * (scalar(1.0) - albGlobe)
        * albSfc * SR;

    scalar qrSun = qrSunDiff + qrSunDir + qrSunRef;

    // residual from energy balance equation
    scalar residual = - (qrG + qhG) + (qrSun + qrAtm + qrSfc);

    scalar rms = sqrt(pow(residual,2));

    return rms;
}


Foam::scalar Foam::comfortFunctions::WBGTiso7243
(
    const scalar& pressure,
    const scalar& va,
    const scalar& Tk,
    const scalar& relativeHumidity,
    const scalar& wetBulbTemperature,
    const scalar& rho,
    const scalar& nu,
    const scalar& Cp,
    const scalar& kappa,
    const scalar& solarRad,
    const Switch& flagSolarRad,
    const scalar& zenithAngle
)
{
    scalar Wbgt = 0;
    scalar vel = va;
    scalar RHfrac = relativeHumidity*0.01;
    scalar Tpwb = wetBulbTemperature;
    scalar SR = solarRad;

    if (va < 0.3)
    {
        vel = 0.3;
    }

    if (flagSolarRad) //outdoor conditions
    {
        // Golden-section Method - find local minimum

        scalar phi = 0.5*(sqrt(5.0) - 1.0);
        scalar lambda = 0.5*(3.0 - sqrt(5.0));

        scalar Tmin = Tk - 4.0;
        scalar Tmax = Tk + 15.0;

            // Tglobe calculation

        scalar xa = Tmin;
        scalar xb = Tmax;
        //scalar fa;
        //scalar fb;

        scalar x1 = xb - phi*(xb - xa);
        scalar x2 = xa + phi*(xb - xa);

        scalar fx1 =
            fTg
            (
                x1, Tk, RHfrac, vel, rho, nu,
                Cp, kappa, SR, zenithAngle
            );

        scalar fx2 =
            fTg
            (
                x2, Tk, RHfrac, vel, rho, nu,
                Cp, kappa, SR, zenithAngle
            );

        while (fabs(xa - xb) > 0.0001)
        {
            if (fx1 > fx2)
            {
                xa = x1;
                //fa = fx1;

                if (fabs(xa - xb) < 0.0001)
                {
                    break;
                }

                x1 = x2;
                fx1 = fx2;
                x2 = xb - lambda*(xb - xa);

                fx2 =
                (
                    fTg
                    (
                        x2, Tk, RHfrac, vel, rho, nu,
                        Cp, kappa, SR, zenithAngle
                    )
                );
            }
            else
            {
                xb = x2;
                //fb = fx2;

                if (fabs(xa - xb) < 0.0001)
                {
                    break;
                }

                x2 = x1;
                fx2 = fx1;
                x1 = xa + lambda*(xb - xa);

                fx1 =
                (
                    fTg
                    (
                        x1, Tk, RHfrac, vel, rho, nu,
                        Cp, kappa, SR, zenithAngle
                    )
                );
            }
        }

        scalar Tglobe = xa;

            // Tnwb calculation

        Tmin = Tk - 20.0;
        Tmax = Tk + 5.0;
        xa = Tmin;
        xb = Tmax;

        x1 = xb - phi*(xb - xa);
        x2 = xa + phi*(xb - xa);

        fx1 =
        (
            fabs
            (
                fTnwb
                (
                    x1, pressure, Tk, RHfrac, vel, rho, nu,
                    Cp, kappa, SR, zenithAngle
                )
              - x1
            )
        );

        fx2 =
        (
            fabs
            (
                fTnwb
                (
                    x2, pressure, Tk, RHfrac, vel, rho, nu,
                    Cp, kappa, SR, zenithAngle
                )
              - x2
            )
        );

        while (fabs(xa - xb) > 0.0001)
        {
            if (fx1 > fx2)
            {
                xa = x1;
                //fa = fx1;

                if (fabs(xa - xb) < 0.0001)
                {
                    break;
                }

                x1 = x2;
                fx1 = fx2;
                x2 = xb - lambda*(xb - xa);

                fx2 =
                (
                    fabs
                    (
                        fTnwb
                        (
                            x2, pressure, Tk, RHfrac, vel, rho, nu,
                            Cp, kappa, SR, zenithAngle
                        )
                      - x2
                    )
                );

            }
            else
            {
                xb = x2;
                //fb = fx2;

                if (fabs(xa - xb) < 0.0001)
                {
                    break;
                }

                x2 = x1;
                fx2 = fx1;
                x1 = xa + lambda*(xb - xa);

                fx1 =
                (
                    fabs
                    (
                        fTnwb
                        (
                            x1, pressure, Tk, RHfrac, vel, rho, nu,
                            Cp, kappa, SR, zenithAngle
                        )
                      - x1
                    )
                );
            }
        }

        scalar Tnwb = xa;

        Wbgt = 0.7*Tnwb + 0.2*Tglobe + 0.1*Tk;

    }
    else // indoor conditions
    {
        if ((vel >= 0.3) && (vel <= 3.0))
        {
            Wbgt = 0.67*Tpwb + 0.33*Tk - 0.048*log(vel)*(Tk - Tpwb);
        }
        else //if (va > 3.0)
        {
            Wbgt = 0.7*Tpwb + 0.3*Tk;
        }
    }

    return Wbgt;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::comfortFunctions::comfortFunctions
(
    const scalar clo,
    const scalar met,
    const scalar work,
    const scalar Tu
)
    :
    clo_(clo),
    met_(met),
    work_(work),
    sigma_(constant::physicoChemical::sigma.value())
{}
