/**
 * Bogus code - never called, but ensures that the OpenBabel formats are referenced.
 *
 * Just needed in Win32, as we're staticly linking rather than dynamic linkage...
 */

#include <openbabel/obmolecformat.h>

/**
 * Note that pointers are all of OBMoleculeFormat - the parent class of each format;
 * this ensures we don't need to drag in each format's definition, which isn't 
 * normally exposed (no header files exist!)
 */
extern OpenBabel::OBMoleculeFormat *theACRFormatPtr;
extern OpenBabel::OBMoleculeFormat *theAlchemyFormatPtr;
extern OpenBabel::OBMoleculeFormat *theAmberPrepFormatPtr;
extern OpenBabel::OBMoleculeFormat *theBallStickFormatPtr;
extern OpenBabel::OBMoleculeFormat *theBGFFormatPtr;
extern OpenBabel::OBMoleculeFormat *theBoxFormatPtr;
extern OpenBabel::OBMoleculeFormat *theCacheFormatPtr;
extern OpenBabel::OBMoleculeFormat *theCANSMIFormatPtr;
extern OpenBabel::OBMoleculeFormat *theMOL2FormatPtr;
extern OpenBabel::OBMoleculeFormat *thePDBFormatPtr;


/**
 * Code is never called - but forced linker & compiler to reference these objects
 * & hence ensures constructor is called by the static object hiding in OpenBabel.
 * This in turn registers the format with OpenBabel (sheesh!)
 */
void touchOpenBabelFormats()
{
    theACRFormatPtr->GetMIMEType();
    theAlchemyFormatPtr->GetMIMEType();
    theAmberPrepFormatPtr->GetMIMEType();
    theBallStickFormatPtr->GetMIMEType();
    theBGFFormatPtr->GetMIMEType();
    theBoxFormatPtr->GetMIMEType();
    theCacheFormatPtr->GetMIMEType();
    theCANSMIFormatPtr->GetMIMEType();
    theMOL2FormatPtr->GetMIMEType();
    thePDBFormatPtr->GetMIMEType();
}
