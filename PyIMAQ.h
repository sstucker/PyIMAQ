// PyIMAQ.h : main header file for the PyIMAQ DLL
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'pch.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CPyIMAQApp
// See PyIMAQ.cpp for the implementation of this class
//

class CPyIMAQApp : public CWinApp
{
public:
	CPyIMAQApp();

// Overrides
public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
};
