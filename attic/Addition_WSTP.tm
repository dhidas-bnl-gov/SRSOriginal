:Begin:
:Function:      OSCARS_Init
:Pattern:       OSCARSInit[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Integer
:End:
:Evaluate: OSCARSInit::usage = "OSCARS_Init[] initialize.  Returns Object number"


:Begin:
:Function:      OSCARS_Delete
:Pattern:       OSCARSDelete[Obj_]
:Arguments:     { Obj }
:ArgumentTypes: { Integer }
:ReturnType:    Integer
:End:
:Evaluate: OSCARSDelete::usage = "OSCARS_Delete[i_] deletes the i_th OSCARS object.  Returns -1 to invalidate reference"


:Begin:
:Function:      OSCARS_DeleteAll
:Pattern:       OSCARSDeleteAll[]
:Arguments:     {  }
:ArgumentTypes: {  }
:ReturnType:    Integer
:End:
:Evaluate: OSCARSDeleteAll::usage = "OSCARS_DeleteAll[] deletes all OSCARS objects.  Returns -1 to indicate success"



:Begin:
:Function:      OSCARS_AddMagneticField
:Pattern:       OSCARSAddMagneticField[Obj_, FileName_, Format_]
:Arguments:     { Obj, FileName, Format }
:ArgumentTypes: { Integer, String, String }
:ReturnType:    Null
:End:
:Evaluate: OSCARSAddMagneticField::usage = "OSCARSAddMagneticField[] adds a magnetic field to OSCARS object Obj_."



:Begin:
:Function:      OSCARS_GetB
:Pattern:       OSCARSGetB[Obj_, Position_List]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Manual
:End:
:Evaluate: OSCARSGetB::usage = "OSCARSGetB[Position_] returns magnetic field at a point"



:Begin:
:Function:      OSCARS_GetBx
:Pattern:       OSCARSGetBx[Obj_, Position_]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Real
:End:
:Evaluate: OSCARSGetBx::usage = "OSCARSGetBx[] returns X-component of magnetic field at a point"



:Begin:
:Function:      OSCARS_GetBy
:Pattern:       OSCARSGetBy[Obj_, Position_]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Real
:End:
:Evaluate: OSCARSGetBy::usage = "OSCARSGetBy[] returns Y-component of magnetic field at a point"



:Begin:
:Function:      OSCARS_GetBz
:Pattern:       OSCARSGetBz[Obj_, Position_]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Real
:End:
:Evaluate: OSCARSGetBz::usage = "OSCARSGetBz[] returns Z-component of magnetic field at a point"






