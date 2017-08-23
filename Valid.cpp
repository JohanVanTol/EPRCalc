//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "Valid.h"
#include <sysutils.hpp>
#include <excpt.h>

//---------------------------------------------------------------------------
#pragma package(smart_init)
// Validation routines
bool ValidReal(AnsiString S, double* pReal)
{
    double d;
    try
    {
        *pReal = S.ToDouble();
    }
    catch (const EConvertError&)
    {
        *pReal = 0;
        return false;
    }
    return true;
}

bool ValidInt(AnsiString S, int* pInt)
{
    try
    {
        *pInt = S.ToDouble();
    }
    catch (const EConvertError&)
    {
        *pInt = 0;
        return false;
    }
    return true;
}
