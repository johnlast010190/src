{
    DDtU1 =
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1;
    fvOptions.addAcceleration(U1, DDtU1);

    DDtU2 =
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2;
    fvOptions.addAcceleration(U2, DDtU2);
}
