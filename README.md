#EULEROFIT
Fortran code to estimate a pole of rotation from GNSS data.
By Nicola D'Agostino, last updated 2024

eulerofit operates in two modes:
- calculates the best-fit Eulerian poles of rotation for the stations defined in file *.fix (format: ID_STATION[4 chars]  Plate[2 Chars]
- subtracts  [rescalc.out] rotation defined in xyzpole.out or enupole.out files
- calculates [velcalc.out] velocities predicted by rotation defined in xyzpole.out or enupole.out files.

The imput data format is the following:

  lon        lat         East        North       sigE        sigN        corrE-N     SiteID
  
 [deg]      [deg]      [mm/yr]      [mm/yr]     [mm/yr]     [mm/yr]      [-1/1]     [4 chars]

eulerofit returns statistical metrics to evaluated the goodness of the fit and plate rigidity.

The fortran code is contained in eulerofit_single_precision.f while eulerofit is a csh driver
to assemble and format input files.

# ------------------------------------------
For usage and an example workflow:

Calculate best fit Eulerian pole of the Apulian microplate using the stations defined in Ap.fix.

./eulerofit -v GMT_Velocita_EU.dat -e Ap.fix


Subtract the Eulerian pole defined in xyzpole.out (Cartesian components).

./eulerofit -v GMT_Velocita_EU.dat -p1 xyzpole.out


Subtract the Eulerian pole defined in enupole.out (lon lat rot).

./eulerofit -v GMT_Velocita_EU.dat -p2 enupole.out


GNSS velocity field from:

Devoti R., N. D'Agostino, E. Serpelloni, G. Pietrantonio, et al., 2017. A Combined Velocity Field of the Mediterranean Region, 

ANNALS OF GEOPHYSICS, 60, 2, 2017, S0217; doi:10.4401/ag-7059.

# ------------------------------------------
To compile eulerofit_single_precision.f:

make -f make_Osx
