    //solution controls
    label msi =
        solverControls.lookupOrDefault<label>("momentumSolutionInterval", 1);

    label tsi =
        solverControls.lookupOrDefault<label>("temperatureSolutionInterval", 1);

    if (solverControls.found("buoyancyFormulation"))
    {
        word buoyancyFormulation(solverControls.lookup("buoyancyFormulation"));

        bool newCellBuoyancy = cellBuoyancy;

        if (buoyancyFormulation == "face")
        {
            newCellBuoyancy = false;
        }
        else if (buoyancyFormulation == "cell")
        {
            newCellBuoyancy = true;
        }
        else
        {
            FatalError << "Unknown buoyancy formulation " << buoyancyFormulation
                << " specified. Valid types are 'cell' and 'face'" << nl
                << exit(FatalError);
        }

        cellBuoyancy = newCellBuoyancy;
    }

    bool flow = true;
    bool thermal = true;

    if (msi > 0)
    {
        flow = (runTime.timeIndex()%msi == 0);
    }
    else
    {
        flow = false;
    }

    if (tsi > 0)
    {
        thermal = (runTime.timeIndex()%tsi == 0);
    }
    else
    {
        thermal = false;
    }
