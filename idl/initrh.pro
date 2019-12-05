;; --- Precompile these analysis routines --              --------- ;

forward_function  openJ, readGeometry, readAtmos, readAtom, $
 openOpacity, getTau, xann, yann, imu, readMolecules, setRatio, $
 readmollines, readInput, timestamp, linecooling, setQuantity

.r misc_func xyann readinput timestamp linecooling
.r readatmos readatom readspectrum readgeometry readopacity
.r getline  readflux  readray  readflow  readmolecules
.r viewatmos viewatom viewpops viewie viewtrans viewmolecules
.r viewmolpops viewsplitting viewstokes readmollines
.r viewtermdiag viewmolterm viewangles viewabundance viewbopac
.r orient viewj viewsource viewavg viewgrid viewdisk viewatlas
.r KPNOatlas Hawaiiatlas KPKatlas sumeratlas KPspotatlas
.r KPIRatlas ATMOSatlas SolFluxatlas KPspotIRatlas
.r rawatom viewalpha viewcontrib xyann imu viewflux rotbroad 
.r rayinterpolate raytrace viewray
.r analyze 

answer = 'n'  &  read, answer, PROMPT="Save RH analysis routines [n]? "
IF (strmid(answer, 0, 1) EQ 'y') THEN $
 save, /ROUTINES, FILENAME=getenv('RH_IDL_PATH')+'/rh.sav'
