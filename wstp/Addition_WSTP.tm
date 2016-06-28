:Begin:
:Function:      SRS_Init
:Pattern:       SRSInit[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Integer
:End:
:Evaluate: SRSInit::usage = "SRS_Init[] initialize.  Returns Object number"


:Begin:
:Function:      SRS_Delete
:Pattern:       SRSDelete[Obj_]
:Arguments:     { Obj }
:ArgumentTypes: { Integer }
:ReturnType:    Integer
:End:
:Evaluate: SRSDelete::usage = "SRS_Delete[i_] deletes the i_th SRS object.  Returns -1 to invalidate reference"


:Begin:
:Function:      SRS_DeleteAll
:Pattern:       SRSDeleteAll[]
:Arguments:     {  }
:ArgumentTypes: {  }
:ReturnType:    Integer
:End:
:Evaluate: SRSDeleteAll::usage = "SRS_DeleteAll[] deletes all SRS objects.  Returns -1 to indicate success"



:Begin:
:Function:      SRS_AddMagneticField
:Pattern:       SRSAddMagneticField[Obj_, FileName_, Format_]
:Arguments:     { Obj, FileName, Format }
:ArgumentTypes: { Integer, String, String }
:ReturnType:    Null
:End:
:Evaluate: SRSAddMagneticField::usage = "SRSAddMagneticField[] adds a magnetic field to SRS object Obj_."



:Begin:
:Function:      SRS_GetB
:Pattern:       SRSGetB[Obj_, Position_List]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Manual
:End:
:Evaluate: SRSGetB::usage = "SRSGetB[Position_] returns magnetic field at a point"



:Begin:
:Function:      SRS_GetBx
:Pattern:       SRSGetBx[Obj_, Position_]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Real
:End:
:Evaluate: SRSGetBx::usage = "SRSGetBx[] returns X-component of magnetic field at a point"



:Begin:
:Function:      SRS_GetBy
:Pattern:       SRSGetBy[Obj_, Position_]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Real
:End:
:Evaluate: SRSGetBy::usage = "SRSGetBy[] returns Y-component of magnetic field at a point"



:Begin:
:Function:      SRS_GetBz
:Pattern:       SRSGetBz[Obj_, Position_]
:Arguments:     { Obj, Position }
:ArgumentTypes: { Integer, RealList }
:ReturnType:    Real
:End:
:Evaluate: SRSGetBz::usage = "SRSGetBz[] returns Z-component of magnetic field at a point"






