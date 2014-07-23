// program_gcm_ds
// Consistent Statistical Downscaling
// for IPCC Global Circulation Models
// Feng Zhang <fzhang@gatech.edu>
// July, 2009

//Revision: 
//July 18, 2009 -- adding downscaling fun.
//July 20, 2009 -- adding n-day aggregation.
//August 10, 2009 -- revising job file structure.
//October 30, 2009 -- adding testing for 20C3M senarios.

#include "stdafx.h"
#include "netcdf_test.h"
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <nag.h>
#include <nag_stdlib.h>
#include <nagg01.h>
#include <nagm01.h>
#include "gwrinc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

//debug options
//#define DEBUG_MODE			  
//#define DEBUG_20C3M
#define DO_BC			      
#define DO_DS			      
#define EXPORT_BC_GRID	      
#define EXPORT_DS_GRID	      
#define OUTPUT_DIAGNOSTIC       

//daily data set constant
#define DO_DAILY_DS           
#define UM_DAILY_SLICEPERMO   6   //30=01D//6=05D//3=10D//1=30D
#define DS_CHECK_BOUND        0   //0-do not ignore;1-ignore out-of-limit DS values

//constant definition
#define MAX_MONTH_NUM         100*12
#define MAX_DAILY_NUM         100*12*UM_DAILY_SLICEPERMO
#define MAX_CHAR_STRLINE      1024
#define MAX_CHAR_PERLINE      1024
#define MAX_WORKING_CELLS     128

//bias correction macros
#define OPT_MARGINE_BC	      0   //0-direct convert;1-ECDF mapping
#define OPT_JOINT_BC	      0   //0-NO 2DECDF maping;1-2DECDF joint maping
#define DATA_COL_NUM          6
#define NUM_BIN               50
#define LMT_JTMG_CDFDIFF      0.1

#ifdef  DEBUG_20C3M
#define DO_ANOMANY            //defined if working on anomany sequences
							  //not defined if working on differentiated sequences
#endif

// The one and only application object
CWinApp theApp;
using namespace std;

//include function header files
#include "gwri_cbc.h"
#include "gwri_cds.h"

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	// initialize MFC and print and error on failure
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: change error code to suit your needs
		_tprintf(_T("Fatal Error: MFC initialization failed\n"));
		nRetCode = 1;
	}
	else
	{
		//------------------------------------------
        //---------- Software Information ----------
        //------------------------------------------
		printf("------------------------------------\n");
		printf("Program Information\n");
		printf("------------------------------------\n");
		printf("Consistent Statistical Downscaling\n");
		printf("<IPCC Global Circulation Models>\n");
		printf("Version: v1.0b, August, 2009\n");
		printf("By Feng Zhang <fzhang@gatech.edu>\n");
		printf("Copyright 2009 Georgia Tech, Atlanta\n");
		printf("\n");

#ifndef DEBUG_MODE
		if(argc != 2)
		{
			printf("Please input the job file name!\n");
			exit(-1);
		}
		CString strJobFilename = argv[1];
#else

		//job files for 20C3M scanarios
		//CCCMA_CGCM31
		//CString strJobFilename = CString("c:\\temp\\pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.20cm.run1.05D.c");

		//job files for SRES_A2 scanarios
		//CCCMA_CGCM31
		//CString strJobFilename = CString("c:\\temp\\pcmdi.ipcc4.cccma_cgcm3_1.sresa2.run1.05D.c");
		//MPI_ECHAM
		CString strJobFilename = CString("c:\\temp\\pcmdi.ipcc4.mpi_echam5.sresa2.run1.05D.c");

		//job files for SRES_A1B scanarios
		//CCCMA_CGCM31
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.05D.c");
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.10D.c");
		
		
		//GFDL_CM21
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.gfdl_cm2_1.sresa1b.run1.05D.c");
		//MPI_ECHAM
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.mpi_echam5.sresa1b.run1.05D.c");
		//MRI_CGCM232
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.mri_cgcm2_3_2a.sresa1b.run1.05D.c");
		//NCAR_CCSM30
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.ncar_ccsm3_0.sresa1b.run2.05D.c");
		//NCAR_PCM1
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.ncar_pcm1.sresa1b.run1.05D.c");
		//UKMO_HADCM3
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.ukmo_hadcm3.sresa1b.run1.05D.c");
		//BCCR_BCM2
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.bccr_bcm2_0.sresa1b.run1.05D.c");
		//CNRM_CM31
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.cnrm_cm3.sresa1b.run1.05D.c");
		//CSIRO_MK3
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.csiro_mk3_0.sresa1b.run1.05D.c");
		//GISS_E_R
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.giss_model_e_r.sresa1b.run2.05D.c");
		//INMCM
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.inmcm3_0.sresa1b.run1.05D.c");
		//MIUB_ECHO_G
		//CString strJobFilename = CString("C:\\temp\\pcmdi.ipcc4.miub_echo_g.sresa1b.run1.05D.c");

#endif
		//------------------------------------------
        //---------- Variable Definition------------
        //------------------------------------------
		//DIAGNOSTIC
		char charInputTSFileName[MAX_CHAR_PERLINE];
		char charOutputTSFileName[MAX_CHAR_PERLINE];
		int  iDiagBCLatIndex=0;
		int  iDiagBCLonIndex=0;
		int  iDiagDSLatIndex=0;
		int  iDiagDSLonIndex=0;

		//INPUT:Observation File Directory(UM)
		char charUMMonthlyDir[MAX_CHAR_PERLINE];
		char charUMDailyDir[MAX_CHAR_PERLINE];

		//INPUT:Downscaling woriking region coordinates
		int iLatLimitMin=0;
		int iLatLimitMax=0;
		int iLonLimitMin=0;
		int iLonLimitMax=0;

		//INPUT:Observation Data Pointers (UM Regrid, Monthly)
		char  *charOBSTavgFileName = NULL;
		char  *charOBSPrcpFileName = NULL;
		float *pOBSLat         = NULL;   int nOBSLatDim;
		float *pOBSLon         = NULL;   int nOBSLonDim;
		int   *pOBSTime        = NULL;   //int nOBSTimeDim;
		float *pOBSData        = NULL;
		char  *charOBSLatName  = "lat";
		char  *charOBSLonName  = "lon";
		int   iOBSStartYear, iOBSEndYear;
		int   iOBSWorkStartYear, iOBSWorkEndYear;
		char  *charOBSVartName = NULL;
		float fOBSMissingValue = 1.e+020f;
		int   iOBSMonthOffset  = 0;

		//INPUT:Observation Data Pointers (UM Regrid, Daily)
		char  *charOBS10DTavgFileName = NULL;
		char  *charOBS10DPrcpFileName = NULL;
		float *pOBS10DLat         = NULL;   int nOBS10DLatDim;
		float *pOBS10DLon         = NULL;   int nOBS10DLonDim;
		int   *pOBS10DTime        = NULL;   //int nOBS10DTimeDim;
		float *pOBS10DData        = NULL;
		char  *charOBS10DLatName  = "lat";
		char  *charOBS10DLonName  = "lon";
		int   iOBS10DStartYear, iOBS10DEndYear;
		int   iOBS10DWorkStartYear, iOBS10DWorkEndYear;
		char  *charOBS10DVartName = NULL;
		float fOBS10DMissingValue = 1.e+020f;
		int   iOBS10DMonthOffset  = 0;

		//INPUT:Observation Data Pointers (UM Original)
		char  *charUMTavgFileName = NULL;
		char  *charUMPrcpFileName = NULL;
		float *pUMLat         = NULL;   int nUMLatDim;
		float *pUMLon         = NULL;   int nUMLonDim;
		int   *pUMTime        = NULL;   //int nUMTimeDim;
		float *pUMData        = NULL;
		float *pUMTData= NULL, *pUMPData= NULL;
		char  *charUMLatName  = "latitude";
		char  *charUMLonName  = "longitude";
		int   iUMStartYear, iUMEndYear;
		int   iUMWorkStartYear, iUMWorkEndYear;
		char  *charUMVartName = NULL;
		float fUMMissingValue = 1.e+020f;
		int   iUMMonthOffset  = 0;

		//INPUT:GCM Data Pointers (Control Run)
		char  *charGCMCONTavgFileName = NULL;
		char  *charGCMCONPrcpFileName = NULL;
		float *pGCMCONLat         = NULL; int nGCMCONLatDim;
		float *pGCMCONLon         = NULL; int nGCMCONLonDim;
		int   *pGCMCONTime        = NULL; //int nGCMCONTimeDim;
		float *pGCMCONData        = NULL;
		char  *charGCMCONLatName  = "lat";
		char  *charGCMCONLonName  = "lon";
		int   iGCMCONStartYear, iGCMCONEndYear;
		int   iGCMCONWorkStartYear, iGCMCONWorkEndYear;
		char  *charGCMCONVartName = NULL;
		float fGCMCONMissingValue = 1.e+020f;
		int   iGCMCONMonthOffset  = 0;

		//INPUT:GCM Data Pointers (Future Run)
		char  *charGCMFUTTavgFileName = NULL;
		char  *charGCMFUTPrcpFileName = NULL;
		float *pGCMFUTLat         = NULL; int nGCMFUTLatDim;
		float *pGCMFUTLon         = NULL; int nGCMFUTLonDim;
		int   *pGCMFUTTime        = NULL; //int nGCMFUTTimeDim;
		float *pGCMFUTData        = NULL;
		char  *charGCMFUTLatName  = "lat";
		char  *charGCMFUTLonName  = "lon";
		int   iGCMFUTStartYear, iGCMFUTEndYear;
		int   iGCMFUTWorkStartYear, iGCMFUTWorkEndYear;
		char  *charGCMFUTVartName = NULL;
		float fGCMFUTMissingValue = 1.e+020f;
		int   iGCMFUTMonthOffset  = 0;

		//OUTPUT:GCM Data Pointer (BC for FUT)
		char  *charGCMBCTavgFileName = NULL;
		char  *charGCMBCPrcpFileName = NULL;
		float *pGCMBCLat         = NULL; int nGCMBCLatDim;
		float *pGCMBCLon         = NULL; int nGCMBCLonDim;
		int   *pGCMBCTime        = NULL; int nGCMBCTimeDim;
		float *pGCMBCData        = NULL;
		char  *charGCMBCLatName  = "lat";
		char  *charGCMBCLonName  = "lon";
		char  *charGCMBCVartName = NULL;
		float fGCMBCMissingValue = 1.e+020f;
		int   iGCMBCMonthOffset  = 0;

		//OUTPUT:GCM Data Pointer (DS for FUT)
		char  *charGCMDSTavgFileName = NULL;
		char  *charGCMDSPrcpFileName = NULL;
		float *pGCMDSLat         = NULL; int nGCMDSLatDim;
		float *pGCMDSLon         = NULL; int nGCMDSLonDim;
		int   *pGCMDSTime        = NULL; int nGCMDSTimeDim;
		float *pGCMDSData        = NULL;
		char  *charDSMBCLatName  = "lat";
		char  *charDSMBCLonName  = "lon";
		char  *charDSMBCVartName = NULL;
		float fGCMDSMissingValue = 1.e+020f;
		int   iGCMDSMonthOffset  = 0;

		//OUTPUT:Data Array
		float* (grid_data_bc_Tavg[MAX_MONTH_NUM]);  int num_month_bc_Tavg = 0;
		float* (grid_data_bc_Prcp[MAX_MONTH_NUM]);  int num_month_bc_Prcp = 0;
		float* (grid_data_dsindex_Tavg[MAX_DAILY_NUM]);
		float* (grid_data_dsindex_Prcp[MAX_DAILY_NUM]);
		float* (grid_data_ds_Tavg[MAX_MONTH_NUM]);  int num_month_ds_Tavg = 0;
		float* (grid_data_ds_Prcp[MAX_MONTH_NUM]);  int num_month_ds_Prcp = 0;

        //------------------------------------------
        //-----------Variable Initialization--------
        //------------------------------------------
		charOBSTavgFileName    = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charOBSPrcpFileName    = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charOBS10DTavgFileName = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charOBS10DPrcpFileName = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMCONTavgFileName = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMCONPrcpFileName = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMFUTTavgFileName = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMFUTPrcpFileName = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMBCTavgFileName  = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMBCPrcpFileName  = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMDSTavgFileName  = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charGCMDSPrcpFileName  = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charUMTavgFileName     = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));
		charUMPrcpFileName     = (char*)malloc(MAX_CHAR_STRLINE*sizeof(char));

        //------------------------------------------
        //Parameters needed to be assigned
		//------------------------------------------
		FILE *pf_input = NULL;
		char charInputParamFileName[MAX_CHAR_STRLINE];
		char charInputLineBuffer[MAX_CHAR_STRLINE];
		strcpy(charInputParamFileName, (LPSTR)(LPCTSTR) strJobFilename);
		if( (pf_input = fopen(charInputParamFileName,"rt")) == NULL)
		{
			printf("Can not open input file %s!\n", charInputParamFileName);
			return -1;
		}

		//Reading Inputs
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		//OBS
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charOBSTavgFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charOBSPrcpFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charOBS10DTavgFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charOBS10DPrcpFileName);
		//CON
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMCONTavgFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMCONPrcpFileName);
		//FUT
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMFUTTavgFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMFUTPrcpFileName);
		//BC
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMBCTavgFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMBCPrcpFileName);
		//DS
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMDSTavgFileName);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charGCMDSPrcpFileName);
		//UM Directory
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charUMMonthlyDir);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charUMDailyDir);
		//Working Date
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%d,%d,%d,%d\n", &iOBSStartYear, &iOBSEndYear, &iOBSWorkStartYear, &iOBSWorkEndYear);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%d,%d,%d,%d\n", &iGCMCONStartYear, &iGCMCONEndYear, &iGCMCONWorkStartYear, &iGCMCONWorkEndYear);
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%d,%d,%d,%d\n", &iGCMFUTStartYear, &iGCMFUTEndYear, &iGCMFUTWorkStartYear, &iGCMFUTWorkEndYear);
		iOBS10DStartYear     = iOBSStartYear;     iOBS10DEndYear     = iOBSEndYear;	
		iOBS10DWorkStartYear = iOBSWorkStartYear; iOBS10DWorkEndYear = iOBSWorkEndYear;
		iUMStartYear     = iOBSStartYear;     iUMEndYear     = iOBSEndYear;	
		iUMWorkStartYear = iOBSWorkStartYear; iUMWorkEndYear = iOBSWorkEndYear;
		//Downscaling Region Information
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%d,%d,%d,%d\n", &iLatLimitMin, &iLatLimitMax, &iLonLimitMin, &iLonLimitMax);
		//Diagnostic Outputfile
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%s\n", charInputTSFileName);
		//Diagnostic output cells' index
		fgets(charInputLineBuffer, MAX_CHAR_STRLINE, pf_input);
		fscanf(pf_input, "%d,%d,%d,%d\n", &iDiagBCLatIndex, &iDiagBCLonIndex, &iDiagDSLatIndex, &iDiagDSLonIndex);

		fclose(pf_input);

		//------------------------------------------
        //---------- Inputs Information ----------
        //------------------------------------------
		printf("----------------------------------\n");
		printf("Job Description\n");
		printf("----------------------------------\n");
		printf("[List of input files]\n");
		printf("OBS Data (Monthly) Directory:\n%s\n",charUMMonthlyDir);
		printf("OBS Data (Daily) Directory:\n%s\n",charUMDailyDir);
		printf("GCM Control Run File(Tavg): \n%s\n",charGCMCONTavgFileName);
		printf("GCM Control Run File(Prcp): \n%s\n",charGCMCONPrcpFileName);
		printf("GCM Projection Run File(Tavg): \n%s\n",charGCMFUTTavgFileName);
		printf("GCM Projection Run File(Prcp): \n%s\n",charGCMFUTPrcpFileName);
		
		printf("\n[List of output files]\n");
		printf("GCM Bias Correcting File(Tavg): \n%s\n",charGCMBCTavgFileName);
		printf("GCM Bias Correcting File(Prcp): \n%s\n",charGCMBCPrcpFileName);
		printf("GCM Downscaling File(Tavg): \n%s\n",charGCMDSTavgFileName);
		printf("GCM Downscaling File(Prcp): \n%s\n",charGCMDSPrcpFileName);
		printf("\n");

		printf("----------------------------------\n");
		printf("Procedures Overview\n");
		printf("----------------------------------\n");
		printf("Step1: Reading Monthly OBS Dataset\n");
		printf("Step2: Reading Dialy OBS Dataset\n");
		printf("Step3: Reading GCM Control Run Dataset\n");
		printf("Step4: Reading GCM Projection Run Dataset\n");
		printf("Step5: Bias Correcting for GCM Projection\n");
		printf("Step6: Downscaling for GCM Projection\n");
		printf("Step7: Writing Bias Correcting Dataset\n");
		printf("Step8: Writing Downscaling Dataset\n");
		printf("\n");
		printf("\n");

#ifdef DEBUG_MODE
		printf("Press 'Enter' to start...\n");
		//getchar();
#endif

		//allocate memory for ts
		float* (fOBSTSTavg[MAX_WORKING_CELLS]); int nOBSTSTavg;
		float* (fOBSTSPrcp[MAX_WORKING_CELLS]); int nOBSTSPrcp;
		float* (fOBS10DTSTavg[MAX_WORKING_CELLS]); int nOBS10DTSTavg;
		float* (fOBS10DTSPrcp[MAX_WORKING_CELLS]); int nOBS10DTSPrcp;
		float* (fCONTSTavg[MAX_WORKING_CELLS]); int nCONTSTavg;
		float* (fCONTSPrcp[MAX_WORKING_CELLS]); int nCONTSPrcp;
		float* (fFUTTSTavg[MAX_WORKING_CELLS]); int nFUTTSTavg;
		float* (fFUTTSPrcp[MAX_WORKING_CELLS]); int nFUTTSPrcp;

		//------------------------------------------
		//--------Read OBS and GCM Grid Data--------
        //------------------------------------------
		
		//Read Upscaled OBS Grid Coordinates
		//-----coordinates-----
		int   nOBSWorkMonth   = 12*(iOBSWorkEndYear-iOBSWorkStartYear+1);
		int   nOBSOffsetMonth = 12*(iOBSWorkStartYear-iOBSStartYear);
		get_grid_coord_dim(charOBSTavgFileName, charOBSLatName, charOBSLonName, &nOBSLatDim, &nOBSLonDim);
		pOBSLat  = (float*) malloc(sizeof(float) * nOBSLatDim);
		pOBSLon  = (float*) malloc(sizeof(float) * nOBSLonDim);
		pOBSData = (float*) malloc(sizeof(float) * nOBSLatDim*nOBSLonDim);
		load_grid_coord(charOBSTavgFileName, charOBSLatName, charOBSLonName, pOBSLat, pOBSLon);
		nOBSTSTavg = nOBSWorkMonth;
		nOBSTSPrcp = nOBSWorkMonth;
		//----calculate recording cells-----
		int iRecordingCellIndex[MAX_WORKING_CELLS];
		int nRecordingCells=0;
		for(int iLatIndex = 0; iLatIndex < nOBSLatDim; iLatIndex++)
		{
			for(int iLonIndex = 0; iLonIndex < nOBSLonDim; iLonIndex++)
			{
				int  iCellIndex = iLatIndex*nOBSLonDim + iLonIndex;
				bool bValueMissing = false;

				//filter out out-of-limit cells
				float fWorkingLat = *(pOBSLat + iLatIndex);
				float fWorkingLon = *(pOBSLon + iLonIndex);
				if(fWorkingLat > (float) iLatLimitMin && fWorkingLat < (float) iLatLimitMax &&
				   fWorkingLon > (float) iLonLimitMin && fWorkingLon < (float) iLonLimitMax)
				{
					iRecordingCellIndex[nRecordingCells] = iCellIndex;
					fOBSTSTavg[nRecordingCells] = (float*)malloc(nOBSTSTavg*sizeof(float));
					fOBSTSPrcp[nRecordingCells] = (float*)malloc(nOBSTSPrcp*sizeof(float));
					nRecordingCells++;
				}
			}
		}

		//Read Upscaled OBS Data (Monthly)
		float *fGridData_obs_Tavg = (float*)malloc(nOBSLatDim*nOBSLonDim*sizeof(float));
		float *fGridData_obs_Prcp = (float*)malloc(nOBSLatDim*nOBSLonDim*sizeof(float));
		printf("Step1: Reading OBS Monthly Dataset\n");
		for(int iMonth = 0; iMonth < nOBSWorkMonth; iMonth++)
		{
			if(iMonth % 50 == 0)
				printf("Reading: %d Months Left.\n", nOBSWorkMonth - iMonth);

			load_netcdf_2dvar(charOBSTavgFileName, "Tavg", iMonth + nOBSOffsetMonth, 
				nOBSLatDim, nOBSLonDim, fGridData_obs_Tavg);
			load_netcdf_2dvar(charOBSPrcpFileName, "Prcp", iMonth + nOBSOffsetMonth, 
				nOBSLatDim, nOBSLonDim, fGridData_obs_Prcp);

			//record obs time sequences monthly
			for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
			{
				float *fpOBSTavgData = fOBSTSTavg[iRecCell];
				float *fpOBSPrcpData = fOBSTSPrcp[iRecCell];
				int    iCellIndex=iRecordingCellIndex[iRecCell];

				*(fpOBSTavgData+iMonth) = *(fGridData_obs_Tavg+iCellIndex);
				*(fpOBSPrcpData+iMonth) = *(fGridData_obs_Prcp+iCellIndex);
			}
		}
		if(!fGridData_obs_Tavg) free(fGridData_obs_Tavg);
		if(!fGridData_obs_Prcp) free(fGridData_obs_Prcp);

		printf("\n");

#ifdef DO_DAILY_DS
		//Read Upscaled OBS Grid Coordinates (Daily)
		int   nOBS10DWorkMonth   = 12*UM_DAILY_SLICEPERMO*(iOBSWorkEndYear-iOBSWorkStartYear+1);
		int   nOBS10DOffsetMonth = 12*UM_DAILY_SLICEPERMO*(iOBSWorkStartYear-iOBSStartYear);
		get_grid_coord_dim(charOBS10DTavgFileName, charOBS10DLatName, charOBS10DLonName, &nOBS10DLatDim, &nOBS10DLonDim);
		pOBS10DLat  = (float*) malloc(sizeof(float) * nOBS10DLatDim);
		pOBS10DLon  = (float*) malloc(sizeof(float) * nOBS10DLonDim);
		pOBS10DData = (float*) malloc(sizeof(float) * nOBS10DLatDim*nOBS10DLonDim);
		load_grid_coord(charOBS10DTavgFileName, charOBS10DLatName, charOBS10DLonName, pOBS10DLat, pOBS10DLon);
		nOBS10DTSTavg = nOBS10DWorkMonth;
		nOBS10DTSPrcp = nOBS10DWorkMonth;
		//check the data compability 
		if(nOBSLatDim!=nOBS10DLatDim)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (lat) are not identical.\n"));
			return -1;
		}
		if(nOBSLonDim!=nOBS10DLonDim)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (lon) are not identical.\n"));
			return -1;
		}

		//Read Upscaled OBS Data (Daily)
		for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
		{
			fOBS10DTSTavg[iRecCell] = (float*)malloc(nOBS10DTSTavg*sizeof(float));
			fOBS10DTSPrcp[iRecCell] = (float*)malloc(nOBS10DTSPrcp*sizeof(float));
		}
		float *fGridData_obsDay_Tavg = (float*)malloc(nOBS10DLatDim*nOBS10DLonDim*sizeof(float));
		float *fGridData_obsDay_Prcp = (float*)malloc(nOBS10DLatDim*nOBS10DLonDim*sizeof(float));
		printf("Step2: Reading OBS Daily Dataset\n");
		for(int iMonth = 0; iMonth < nOBS10DWorkMonth; iMonth++)
		{
			if(iMonth % 50 == 0)
				printf("Reading: %d Months Left.\n", nOBS10DWorkMonth - iMonth);
			
			load_netcdf_2dvar(charOBS10DTavgFileName, "Tavg", iMonth + nOBS10DOffsetMonth, 
				nOBS10DLatDim, nOBS10DLonDim, fGridData_obsDay_Tavg);
			load_netcdf_2dvar(charOBS10DPrcpFileName, "Prcp", iMonth + nOBS10DOffsetMonth, 
				nOBS10DLatDim, nOBS10DLonDim, fGridData_obsDay_Prcp);

			//record obs_days time sequences monthly
			for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
			{
				float *fpOBSTavgData = fOBS10DTSTavg[iRecCell];
				float *fpOBSPrcpData = fOBS10DTSPrcp[iRecCell];
				int    iCellIndex=iRecordingCellIndex[iRecCell];

				*(fpOBSTavgData+iMonth) = *(fGridData_obsDay_Tavg+iCellIndex);
				*(fpOBSPrcpData+iMonth) = *(fGridData_obsDay_Prcp+iCellIndex);
			}
		}
		if(!fGridData_obsDay_Tavg) free(fGridData_obsDay_Tavg);
		if(!fGridData_obsDay_Prcp) free(fGridData_obsDay_Prcp);
		printf("\n");
#endif

		//Read Original OBS Grid Coordinates
		sprintf(charUMTavgFileName, charUMMonthlyDir, "Tavg", iOBS10DStartYear);
		int   nUMWorkMonth   = 12*(iUMWorkEndYear-iUMWorkStartYear+1);
		int   nUMOffsetMonth = 12*(iUMWorkStartYear-iUMStartYear);
		get_grid_coord_dim(charUMTavgFileName, charUMLatName, charUMLonName, &nUMLatDim, &nUMLonDim);
		pUMLat  = (float*) malloc(sizeof(float) * nUMLatDim);
		pUMLon  = (float*) malloc(sizeof(float) * nUMLonDim);
		pUMData = (float*) malloc(sizeof(float) * nUMLatDim*nUMLonDim);
		pUMTData = (float*) malloc(sizeof(float) * nUMLatDim*nUMLonDim);
		pUMPData = (float*) malloc(sizeof(float) * nUMLatDim*nUMLonDim);
		load_grid_coord(charUMTavgFileName, charUMLatName, charUMLonName, pUMLat, pUMLon);
		printf("\n");

		//Read GCM_CON Grid Coordinates
		int   nGCMCONWorkMonth   = 12*(iGCMCONWorkEndYear-iGCMCONWorkStartYear+1);
		int   nGCMCONOffsetMonth = 12*(iGCMCONWorkStartYear-iGCMCONStartYear);
		get_grid_coord_dim(charGCMCONTavgFileName, charGCMCONLatName, charGCMCONLonName, &nGCMCONLatDim, &nGCMCONLonDim);
		pGCMCONLat  = (float*) malloc(sizeof(float) * nGCMCONLatDim);
		pGCMCONLon  = (float*) malloc(sizeof(float) * nGCMCONLonDim);
		pGCMCONData = (float*) malloc(sizeof(float) * nGCMCONLatDim*nGCMCONLonDim);
		load_grid_coord(charGCMCONTavgFileName, charGCMCONLatName, charGCMCONLonName, pGCMCONLat, pGCMCONLon);
		nCONTSTavg = nGCMCONWorkMonth;
		nCONTSPrcp = nGCMCONWorkMonth;
		//check the data compability 
		if(nOBSLatDim!=nGCMCONLatDim)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (lat) are not identical.\n"));
			return -1;
		}
		if(nOBSLonDim!=nGCMCONLonDim)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (lon) are not identical.\n"));
			return -1;
		}
		if(nOBSWorkMonth!=nGCMCONWorkMonth)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (time) are not identical.\n"));
			return -1;
		}

		//Read GCM_CON Grid Data
		for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
		{
			fCONTSTavg[iRecCell] = (float*)malloc(nCONTSTavg*sizeof(float));
			fCONTSPrcp[iRecCell] = (float*)malloc(nCONTSPrcp*sizeof(float));
		}
		float *fGridData_con_Tavg = (float*)malloc(nGCMCONLatDim*nGCMCONLonDim*sizeof(float));
		float *fGridData_con_Prcp = (float*)malloc(nGCMCONLatDim*nGCMCONLonDim*sizeof(float));
		printf("Step3: Reading GCM Control Run Dataset\n");
		for(int iMonth = 0; iMonth < nGCMCONWorkMonth; iMonth++)
		{
			if(iMonth % 50 == 0)
				printf("Reading: %d Months Left.\n", nGCMCONWorkMonth - iMonth);

			load_netcdf_2dvar(charGCMCONTavgFileName, "tas", iMonth + nGCMCONOffsetMonth, 
				nGCMCONLatDim, nGCMCONLonDim, fGridData_con_Tavg);
			load_netcdf_2dvar(charGCMCONPrcpFileName, "pr", iMonth + nGCMCONOffsetMonth, 
				nGCMCONLatDim, nGCMCONLonDim, fGridData_con_Prcp);
			
			//record con time sequences monthly
			for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
			{
				float *fpCONTavgData = fCONTSTavg[iRecCell];
				float *fpCONPrcpData = fCONTSPrcp[iRecCell];
				int    iCellIndex=iRecordingCellIndex[iRecCell];

				*(fpCONTavgData+iMonth) = *(fGridData_con_Tavg+iCellIndex);
				*(fpCONPrcpData+iMonth) = *(fGridData_con_Prcp+iCellIndex);
			}
		}
		if(!fGridData_con_Tavg) free(fGridData_con_Tavg);
		if(!fGridData_con_Prcp) free(fGridData_con_Prcp);
		printf("\n");

		//Read GCM_FUT Grid Coordinates
		int   nGCMFUTWorkMonth   = 12*(iGCMFUTWorkEndYear-iGCMFUTWorkStartYear+1);
		int   nGCMFUTOffsetMonth = 12*(iGCMFUTWorkStartYear-iGCMFUTStartYear);
		get_grid_coord_dim(charGCMFUTTavgFileName, charGCMFUTLatName, charGCMFUTLonName, &nGCMFUTLatDim, &nGCMFUTLonDim);
		pGCMFUTLat  = (float*) malloc(sizeof(float) * nGCMFUTLatDim);
		pGCMFUTLon  = (float*) malloc(sizeof(float) * nGCMFUTLonDim);
		pGCMFUTData = (float*) malloc(sizeof(float) * nGCMFUTLatDim*nGCMFUTLonDim);
		load_grid_coord(charGCMFUTTavgFileName, charGCMFUTLatName, charGCMFUTLonName, pGCMFUTLat, pGCMFUTLon);
		nFUTTSTavg = nGCMFUTWorkMonth;
		nFUTTSPrcp = nGCMFUTWorkMonth;
		//check the data compability 
		if(nOBSLatDim!=nGCMFUTLatDim)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (lat) are not identical.\n"));
			return -1;
		}
		if(nOBSLonDim!=nGCMFUTLonDim)
		{
			_tprintf(_T("Fatal Error: Grid dimensions (lon) are not identical.\n"));
			return -1;
		}
		
		//Read GCM_FUT Grid Data
		//Initialize Data Array of BC+DS
		for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
		{
			fFUTTSTavg[iRecCell] = (float*)malloc(nFUTTSTavg*sizeof(float));
			fFUTTSPrcp[iRecCell] = (float*)malloc(nFUTTSPrcp*sizeof(float));
		}
		float *fGridData_fut_Tavg = (float*)malloc(nGCMFUTLatDim*nGCMFUTLonDim*sizeof(float));
		float *fGridData_fut_Prcp = (float*)malloc(nGCMFUTLatDim*nGCMFUTLonDim*sizeof(float));
		printf("Step4: Reading GCM Projection Run Dataset\n");
		for(int iMonth = 0; iMonth < nGCMFUTWorkMonth; iMonth++)
		{
			if(iMonth % 50 == 0)
				printf("Reading: %d Months Left.\n", nGCMFUTWorkMonth - iMonth);

			load_netcdf_2dvar(charGCMFUTTavgFileName, "tas", iMonth + nGCMFUTOffsetMonth, 
				nGCMFUTLatDim, nGCMFUTLonDim, fGridData_fut_Tavg);
			load_netcdf_2dvar(charGCMFUTPrcpFileName, "pr", iMonth + nGCMFUTOffsetMonth, 
				nGCMFUTLatDim, nGCMFUTLonDim, fGridData_fut_Prcp);

			//record con time sequences monthly
			for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
			{
				float *fpFUTTavgData = fFUTTSTavg[iRecCell];
				float *fpFUTPrcpData = fFUTTSPrcp[iRecCell];
				int    iCellIndex=iRecordingCellIndex[iRecCell];

				*(fpFUTTavgData+iMonth) = *(fGridData_fut_Tavg+iCellIndex);
				*(fpFUTPrcpData+iMonth) = *(fGridData_fut_Prcp+iCellIndex);
			}

			grid_data_bc_Tavg[iMonth]  = (float*)malloc(nGCMFUTLatDim*nGCMFUTLonDim*sizeof(float));
			grid_data_dsindex_Tavg[iMonth]  = (float*)malloc(nGCMFUTLatDim*nGCMFUTLonDim*sizeof(float));
			grid_data_ds_Tavg[iMonth]  = (float*)malloc(nUMLatDim*nUMLonDim*sizeof(float));
			grid_data_bc_Prcp[iMonth]  = (float*)malloc(nGCMFUTLatDim*nGCMFUTLonDim*sizeof(float));
			grid_data_dsindex_Prcp[iMonth]  = (float*)malloc(nGCMFUTLatDim*nGCMFUTLonDim*sizeof(float));
			grid_data_ds_Prcp[iMonth]  = (float*)malloc(nUMLatDim*nUMLonDim*sizeof(float));
			num_month_bc_Tavg ++;
			num_month_ds_Tavg ++;
			num_month_bc_Prcp ++;
			num_month_ds_Prcp ++;
		}
		if(!fGridData_fut_Tavg) free(fGridData_fut_Tavg);
		if(!fGridData_fut_Prcp) free(fGridData_fut_Prcp);
		printf("\n");

#ifdef DO_BC
		//------------------------------------------
        //--------Bias Correct Cell by Cell---------
        //------------------------------------------
		
		//allocate memory for ts
		float* fOBSTSTavgTS;
		float* fOBSTSPrcpTS;
		float* fOBS10DTSTavgTS;
		float* fOBS10DOutTavgTS;
		float* fOBS10DTSPrcpTS;
		float* fOBS10DOutPrcpTS;
		float* fCONTSTavgTS;
		float* fCONTSPrcpTS;
		float* fFUTTSTavgTS;
		float* fFUTTSPrcpTS;
		fOBSTSTavgTS = (float*)malloc(nOBSTSTavg*sizeof(float));
		fOBSTSPrcpTS = (float*)malloc(nOBSTSPrcp*sizeof(float));
		fOBS10DTSTavgTS = (float*)malloc(nOBS10DTSTavg*sizeof(float));
		fOBS10DOutTavgTS = (float*)malloc(nOBS10DTSTavg*sizeof(float));
		fOBS10DTSPrcpTS = (float*)malloc(nOBS10DTSPrcp*sizeof(float));
		fOBS10DOutPrcpTS = (float*)malloc(nOBS10DTSPrcp*sizeof(float));
		fCONTSTavgTS = (float*)malloc(nCONTSTavg*sizeof(float));
		fCONTSPrcpTS = (float*)malloc(nCONTSPrcp*sizeof(float));
		fFUTTSTavgTS = (float*)malloc(nFUTTSTavg*sizeof(float));
		fFUTTSPrcpTS = (float*)malloc(nFUTTSPrcp*sizeof(float));
		
		float* fBCTSTavgTS;  int nBCTSTavg;
		float* fBCTSPrcpTS;  int nBCTSPrcp;
		float* fIndexDSTSTavgTS;  int nDSTSTavg;
		float* fIndexDSTSPrcpTS;  int nDSTSPrcp;
		nBCTSTavg      = nGCMFUTWorkMonth;
		fBCTSTavgTS      = (float*)malloc(nBCTSTavg*sizeof(float));
		nBCTSPrcp      = nGCMFUTWorkMonth;
		fBCTSPrcpTS      = (float*)malloc(nBCTSPrcp*sizeof(float));
		nDSTSTavg      = nGCMFUTWorkMonth;
		fIndexDSTSTavgTS = (float*)malloc(nBCTSTavg*sizeof(float));
		nDSTSPrcp      = nGCMFUTWorkMonth;
		fIndexDSTSPrcpTS = (float*)malloc(nBCTSPrcp*sizeof(float));

		int iOutLon = iDiagBCLonIndex; int iOutLat = iDiagBCLatIndex;
		int iCellOutIndex = iOutLat*nGCMFUTLonDim + iOutLon;

		printf("Step5: Bias Correcting for GCM Projection\n");
		int iWorkingCellIndex[MAX_WORKING_CELLS];
		int nWorkingCells=0;
		
		//initializing bc array(grid)
		for(int iLatIndex = 0; iLatIndex < nOBSLatDim; iLatIndex++)
		{
			for(int iLonIndex = 0; iLonIndex < nOBSLonDim; iLonIndex++)
			{
				int  iCellIndex = iLatIndex*nOBSLonDim + iLonIndex;
				for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
				{
					float *fpBCTavgData = grid_data_bc_Tavg[iMonth];
					float *fpBCPrcpData = grid_data_bc_Prcp[iMonth];
					*(fpBCTavgData+iCellIndex) = fGCMBCMissingValue;
					*(fpBCPrcpData+iCellIndex) = fGCMBCMissingValue;
				}
			}
		}

		for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
		{
			int  iCellIndex=iRecordingCellIndex[iRecCell];
			bool bValueMissing = false;

			//retrieve obs time sequences monthly
			for(int iMonth = 0; iMonth < nOBSTSTavg; iMonth++)
			{
				float *fpTavgData = fOBSTSTavg[iRecCell];
				float *fpPrcpData = fOBSTSPrcp[iRecCell];
				float fTempTavg = *(fpTavgData+iMonth);
				float fTempPrcp = *(fpPrcpData+iMonth);

				if(fTempTavg == fOBSMissingValue || fTempPrcp == fOBSMissingValue)
				{
					bValueMissing = true;
					break;
				}
				*(fOBSTSTavgTS+iMonth) = fTempTavg;
				*(fOBSTSPrcpTS+iMonth) = fTempPrcp;
			}

#ifdef DO_DAILY_DS
			//retrieve obs time sequences daily
			for(int iMonth = 0; iMonth < nOBS10DTSTavg; iMonth++)
			{
				float *fpTavgData = fOBS10DTSTavg[iRecCell];
				float *fpPrcpData = fOBS10DTSPrcp[iRecCell];
				float fTempTavg = *(fpTavgData+iMonth);
				float fTempPrcp = *(fpPrcpData+iMonth);

				if(fTempTavg == fOBS10DMissingValue || fTempPrcp == fOBS10DMissingValue)
				{
					bValueMissing = true;
					break;
				}
				*(fOBS10DTSTavgTS+iMonth) = fTempTavg;
				*(fOBS10DTSPrcpTS+iMonth) = fTempPrcp;

				if(iCellOutIndex == iCellIndex)
				{
					*(fOBS10DOutTavgTS+iMonth) = fTempTavg;
					*(fOBS10DOutPrcpTS+iMonth) = fTempPrcp;
				}
			}
#endif

			//retrieve gcm con time sequences
			for(int iMonth = 0; iMonth < nCONTSTavg; iMonth++)
			{
				float *fpTavgData = fCONTSTavg[iRecCell];
				float *fpPrcpData = fCONTSPrcp[iRecCell];
				float fTempTavg = *(fpTavgData+iMonth);
				float fTempPrcp = *(fpPrcpData+iMonth);

				if(fTempTavg != fGCMCONMissingValue && fTempTavg < 1e5)
				{
					*(fCONTSTavgTS+iMonth) = (float)(fTempTavg-273.15); //convert K to C;
				}
				else
				{
					*(fCONTSTavgTS+iMonth) = fGCMCONMissingValue;
					bValueMissing = true;
					break;
				}

				if(fTempPrcp != fGCMCONMissingValue && fTempPrcp < 1e5)
				{
					*(fCONTSPrcpTS+iMonth) = fTempPrcp*3600*24; //convert kg/m^2/s to mm/day;
				}
				else
				{
					*(fCONTSPrcpTS+iMonth) = fGCMCONMissingValue;
					bValueMissing = true;
					break;
				}
			}

			//retrieve gcm fut time sequences
			for(int iMonth = 0; iMonth < nFUTTSTavg; iMonth++)
			{
				float *fpTavgData = fFUTTSTavg[iRecCell];
				float *fpPrcpData = fFUTTSPrcp[iRecCell];
				float fTempTavg = *(fpTavgData+iMonth);
				float fTempPrcp = *(fpPrcpData+iMonth);

				if(fTempTavg != fGCMFUTMissingValue && fTempTavg < 1e5)
				{
					*(fFUTTSTavgTS+iMonth) = (float)(fTempTavg-273.15); //Convert K to C;
				}
				else
				{
					*(fFUTTSTavgTS+iMonth) = fGCMFUTMissingValue;
					bValueMissing = true;
					break;
				}

				if(fTempPrcp != fGCMFUTMissingValue && fTempPrcp < 1e5)
				{
					*(fFUTTSPrcpTS+iMonth) = fTempPrcp*3600*24; //Convert kg/m^2/s to mm/day;
				}
				else
				{
					*(fFUTTSPrcpTS+iMonth) = fGCMFUTMissingValue;
					bValueMissing = true;
					break;
				}
			}

			//if Missing values found, do not do bias Correcting
			//but skip this cell
			if(bValueMissing)
			{
				//write Missing value to bc array
				for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
				{
					float *fpBCTavgData = grid_data_bc_Tavg[iMonth];
					float *fpBCPrcpData = grid_data_bc_Prcp[iMonth];
					*(fpBCTavgData+iCellIndex) = fGCMBCMissingValue;
					*(fpBCPrcpData+iCellIndex) = fGCMBCMissingValue;
				}
				continue;
			}

			int   iLatIndex   = int(iCellIndex/nOBSLonDim);
			int   iLonIndex   = iCellIndex%nOBSLonDim;
			float fWorkingLat = *(pOBSLat + iLatIndex);
			float fWorkingLon = *(pOBSLon + iLonIndex);

			printf("Working on Cell %d.\n", iCellIndex);
			printf("Latitute = %f(%d), Longitute = %f(%d).\n",fWorkingLat, iLatIndex, fWorkingLon, iLonIndex);
			iWorkingCellIndex[nWorkingCells] = iCellIndex;
			nWorkingCells++;

			//call ts_bc function to do 
			//bias Correcting for such cell
#ifdef DO_ANOMANY
			int ret = ts_bc_ecdf2(fOBSTSTavgTS, nOBSTSTavg,
					   fOBSTSPrcpTS, nOBSTSPrcp,
					   fCONTSTavgTS, nCONTSTavg,
					   fCONTSPrcpTS, nCONTSPrcp,
					   fFUTTSTavgTS, nFUTTSTavg,
					   fFUTTSPrcpTS, nFUTTSPrcp,
					   fBCTSTavgTS,  nBCTSTavg,
					   fBCTSPrcpTS,  nBCTSPrcp,
					   OPT_MARGINE_BC, OPT_JOINT_BC);
#else
			int ret = ts_bc_ecdf(fOBSTSTavgTS, nOBSTSTavg,
					   fOBSTSPrcpTS, nOBSTSPrcp,
					   fCONTSTavgTS, nCONTSTavg,
					   fCONTSPrcpTS, nCONTSPrcp,
					   fFUTTSTavgTS, nFUTTSTavg,
					   fFUTTSPrcpTS, nFUTTSPrcp,
					   fBCTSTavgTS,  nBCTSTavg,
					   fBCTSPrcpTS,  nBCTSPrcp,
					   OPT_MARGINE_BC, OPT_JOINT_BC);
#endif
			if(ret < 0)
			{
				_tprintf(_T("Fatal Error: ts_bc_ecdf.\n"));
				return -1;
			}

			//write to bc array
			for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
			{
				float *fpBCTavgData = grid_data_bc_Tavg[iMonth];
				float *fpBCPrcpData = grid_data_bc_Prcp[iMonth];
				*(fpBCTavgData+iCellIndex) = *(fBCTSTavgTS+iMonth);
				*(fpBCPrcpData+iCellIndex) = *(fBCTSPrcpTS+iMonth);
			}

			//call ds_search function to do
			//downscaling searching (monthly sequences)
			//and return the index of monthly sequences
			ret = ds_search_monthly(fOBSTSTavgTS,      nOBSTSTavg,
					  fOBSTSPrcpTS,      nOBSTSPrcp,
					  fBCTSTavgTS,       nBCTSTavg,
					  fBCTSPrcpTS,       nBCTSPrcp,
					  fIndexDSTSTavgTS,  nDSTSTavg,
					  fIndexDSTSPrcpTS,  nDSTSPrcp,
					  1, fGCMBCMissingValue);
			if(ret < 0)
			{
				_tprintf(_T("Fatal Error: ds_search_monthly.\n"));
				return -1;
			}

#ifdef DO_DAILY_DS

			//call ds_search function to do
			//downscaling searching (10-daily sequences)
			//and return the index of 10-daily sequences
			//the indexes returned are negtives
			ret = ds_search_daily(fOBS10DTSTavgTS,      nOBS10DTSTavg,
					  fOBS10DTSPrcpTS,      nOBS10DTSPrcp,
					  fBCTSTavgTS,       nBCTSTavg,
					  fBCTSPrcpTS,       nBCTSPrcp,
					  fIndexDSTSTavgTS,  nDSTSTavg,
					  fIndexDSTSPrcpTS,  nDSTSPrcp,
					  1, fGCMBCMissingValue, DS_CHECK_BOUND);
			if(ret < 0)
			{
				_tprintf(_T("Fatal Error: ds_search_daily.\n"));
				return -1;
			}
#endif

			//write to ds_index array
			for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
			{
				float *fpDSITavgData = grid_data_dsindex_Tavg[iMonth];
				float *fpDSIPrcpData = grid_data_dsindex_Prcp[iMonth];
				*(fpDSITavgData+iCellIndex) = *(fIndexDSTSTavgTS+iMonth);
				*(fpDSIPrcpData+iCellIndex) = *(fIndexDSTSPrcpTS+iMonth);
			}


		}//for(iRecCell)

		//free allocated memory
		free(fOBSTSTavgTS);
		free(fOBSTSPrcpTS);
		free(fCONTSTavgTS);
		free(fCONTSPrcpTS);
		free(fFUTTSTavgTS);
		free(fFUTTSPrcpTS);
		free(fBCTSTavgTS);
		free(fBCTSPrcpTS);
		free(fIndexDSTSTavgTS);
		free(fIndexDSTSPrcpTS);
		printf("\n");
		printf("\n");

        //------------------------------------------
        //--------------Reindex ds_index array------
        //------------------------------------------
		if(nWorkingCells==0)
		{
			printf("No working cells are found!\n");
			return 1;
		}
		for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
		{
#ifdef DEBUG_20C3M					
			int iRandomCellIndex = iWorkingCellIndex[0];
#else
			int iRandomCellIndex = iWorkingCellIndex[rand() % nWorkingCells];
#endif
			for(int iWorkCell=0; iWorkCell<nWorkingCells; iWorkCell++)
			{
				int iCellIndex = iWorkingCellIndex[iWorkCell];
				float *fpDSITavgData = grid_data_dsindex_Tavg[iMonth];
				float *fpDSIPrcpData = grid_data_dsindex_Prcp[iMonth];
				*(fpDSITavgData+iCellIndex) = *(fpDSITavgData+iRandomCellIndex);
				*(fpDSIPrcpData+iCellIndex) = *(fpDSIPrcpData+iRandomCellIndex);
			}
		}

#endif

#ifdef DO_DS
        //------------------------------------------
        //--------------Retirve DS Mapping----------
        //------------------------------------------
		printf("Step6: Downscaling for GCM Projection.\nDownscaling:[Total Cells=%d]\n", nWorkingCells);
		//initializing ds array(grid)
		for(int iLatIndex = 0; iLatIndex < nUMLatDim; iLatIndex++)
		{
			for(int iLonIndex = 0; iLonIndex < nUMLonDim; iLonIndex++)
			{
				int  iCellIndex = iLatIndex*nUMLonDim + iLonIndex;
				for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
				{
					float *fpDSTavgData = grid_data_ds_Tavg[iMonth];
					float *fpDSPrcpData = grid_data_ds_Prcp[iMonth];
					*(fpDSTavgData + iCellIndex) = fGCMDSMissingValue;
					*(fpDSPrcpData + iCellIndex) = fGCMDSMissingValue;
				}
			}
		}

		//look over each cell
		int iProcessCount = 0;
		for(int iLatIndex = 0; iLatIndex < nOBSLatDim; iLatIndex++)
		{
			for(int iLonIndex = 0; iLonIndex < nOBSLonDim; iLonIndex++)
			{
				int  iCellIndex = iLatIndex*nOBSLonDim + iLonIndex;
				bool bMissingValueCell = false;

				//ignore cells with 
				//missing values filled in
				for(int iMonth = 0; iMonth < 12; iMonth++)
				{//only check the first year now
					float *fpBCTavgData = grid_data_bc_Tavg[iMonth];
					float *fpBCPrcpData = grid_data_bc_Prcp[iMonth];

					float fBCT = *(fpBCTavgData+iCellIndex);
					float fBCP = *(fpBCPrcpData+iCellIndex);

					if(fBCT == fGCMBCMissingValue || fBCP == fGCMBCMissingValue)
					{
						bMissingValueCell = true;
						continue;
					}
				}
				if(bMissingValueCell)
					continue;

				iProcessCount++;
				printf("Processing Cell#=%d\n", iProcessCount);
				//get lat and lon 
				//on GCM grid
				float fGCMLat = *(pOBSLat+iLatIndex);
				float fGCMLon = *(pOBSLon+iLonIndex);

				//get overlapping cell indexes 
				//on UM OBS grid (high resolution)
				int   npMaxOverlapLat = 0, npMinOverlapLat = 0;
				int   npMaxOverlapLon = 0, npMinOverlapLon = 0;
				get_overlap_cells(pUMLat,  nUMLatDim,
					              pUMLon,  nUMLonDim,
								  pOBSLat, nOBSLatDim,
					              pOBSLon, nOBSLonDim,
								  fGCMLat, fGCMLon,
								  &npMaxOverlapLat, &npMinOverlapLat,
								  &npMaxOverlapLon, &npMinOverlapLon, 0.6f);

				//look over each time slice
				for(int iMonth = 0; iMonth < nBCTSTavg; iMonth++)
				{
#ifndef DEBUG_20C3M					
					if(iMonth % 50 == 0)
						printf("Processing: %d Months Left.\n", nBCTSTavg - iMonth);
#endif

					float *fpDSITavgData = grid_data_dsindex_Tavg[iMonth];
					float *fpDSIPrcpData = grid_data_dsindex_Prcp[iMonth];

					//get ds index
					float fDSIT = *(fpDSITavgData+iCellIndex);
					float fDSIP = *(fpDSIPrcpData+iCellIndex);

					//outof limit cell
					//ignore it
					if(fDSIT == fGCMBCMissingValue || fDSIP ==  fGCMBCMissingValue)
						continue;

					//convert them to int type
					int   iDSIT = (int)fDSIT;
					int   iDSIP = (int)fDSIP;

					//negtive values mean
					//daily sequences
					if(fDSIT < 0 || fDSIP < 0)
					{
						//look into daily observational dataset
						//get daily OBS file name
						iDSIT = -iDSIT;
						iDSIP = -iDSIP;
						int   iDSTYear  = int(iDSIT/(12*UM_DAILY_SLICEPERMO)) + iOBSWorkStartYear;
						int   iDSTMonth = iDSIT%(12*UM_DAILY_SLICEPERMO);
						int   iDSPYear  = int(iDSIP/(12*UM_DAILY_SLICEPERMO)) + iOBSWorkStartYear;
						int   iDSPMonth = iDSIP%(12*UM_DAILY_SLICEPERMO);

						char  charMonthlyOBSTavgFileName[MAX_CHAR_STRLINE];
						char  charMonthlyOBSPrcpFileName[MAX_CHAR_STRLINE];
						sprintf(charMonthlyOBSTavgFileName, charUMDailyDir, "Tavg", iDSTYear);
						sprintf(charMonthlyOBSPrcpFileName, charUMDailyDir, "Prcp", iDSPYear);

						//load one time slice Tavg data
						//load one time slice Prcp data
						float *fpDSTavgData = grid_data_ds_Tavg[iMonth];
						float *fpDSPrcpData = grid_data_ds_Prcp[iMonth];

						int   iCount = 0;
						int   iDSTimeSlice=iDSTMonth*(30/UM_DAILY_SLICEPERMO);
						for(int iTimeSlice=iDSTimeSlice; iTimeSlice<iDSTimeSlice+(30/UM_DAILY_SLICEPERMO); iTimeSlice++)
						{
							//printf("%s\n", charMonthlyOBSTavgFileName);
							load_netcdf_2dvar(charMonthlyOBSTavgFileName, "Tavg", iTimeSlice, nUMLatDim, nUMLonDim, pUMTData);
							load_netcdf_2dvar(charMonthlyOBSPrcpFileName, "Prcp", iTimeSlice, nUMLatDim, nUMLonDim, pUMPData);
						
							//copy data from UM OBS to DS Result Array
							for(int iHighLatIndex = npMinOverlapLat; iHighLatIndex < npMaxOverlapLat; iHighLatIndex++)
							{
								for(int iHighLonIndex = npMinOverlapLon; iHighLonIndex < npMaxOverlapLon; iHighLonIndex++)
								{
									int  iCellIndex2 = iHighLatIndex*nUMLonDim+iHighLonIndex;

									float fRetriveT = *(pUMTData+iCellIndex2);
									float fRetriveP = *(pUMPData+iCellIndex2);

									if(fRetriveT == fGCMBCMissingValue || fRetriveP == fGCMBCMissingValue)
									{
										*(fpDSTavgData+iCellIndex2) = fGCMBCMissingValue;
										*(fpDSPrcpData+iCellIndex2) = fGCMBCMissingValue;
										continue;
									}
									
									if(iCount == 0)
									{
										*(fpDSTavgData+iCellIndex2) = fRetriveT;
										*(fpDSPrcpData+iCellIndex2) = fRetriveP;
									}
									else
									{
										*(fpDSTavgData+iCellIndex2) = *(fpDSTavgData+iCellIndex2) + fRetriveT;
										*(fpDSPrcpData+iCellIndex2) = *(fpDSPrcpData+iCellIndex2) + fRetriveP;
									}

									if(iTimeSlice == iDSTimeSlice+(30/UM_DAILY_SLICEPERMO)-1)
									{
										*(fpDSTavgData+iCellIndex2) = *(fpDSTavgData+iCellIndex2)/(iCount+1);
										*(fpDSPrcpData+iCellIndex2) = *(fpDSPrcpData+iCellIndex2)/(iCount+1);
									}

								}//for(iHighLatIndex)
							}//for(iHighLonIndex)

							iCount++;
						}//for(iTimeSlice)

						continue;
					}

					//get monthly OBS file name 
					//and month id = 0-11
					int  iDSTYear  = int(iDSIT/12) + iOBSWorkStartYear;
					int  iDSTMonth = iDSIT%12;
					int  iDSPYear  = int(iDSIP/12) + iOBSWorkStartYear;
					int  iDSPMonth = iDSIP%12;
					char charMonthlyOBSTavgFileName[MAX_CHAR_STRLINE];
					char charMonthlyOBSPrcpFileName[MAX_CHAR_STRLINE];

/////////////////////////////////////////	
#ifdef DEBUG_20C3M					
					if(iMonth%12 == 0)
					{
						printf("Working on Month=%02d\t", (iMonth%12)+1);
						printf("Year=%d, Month=%02d.\n",iDSPYear,iDSPMonth+1);
					}
#endif
/////////////////////////////////////////					


					sprintf(charMonthlyOBSTavgFileName, charUMMonthlyDir, "Tavg", iDSTYear);
					sprintf(charMonthlyOBSPrcpFileName, charUMMonthlyDir, "Prcp", iDSPYear);

					//load one time slice Tavg data
					//load one time slice Prcp data
					float *fpDSTavgData = grid_data_ds_Tavg[iMonth];
					float *fpDSPrcpData = grid_data_ds_Prcp[iMonth];
					//printf("%s\n", charMonthlyOBSTavgFileName);
					load_netcdf_2dvar(charMonthlyOBSTavgFileName, "Tavg", iDSTMonth, nUMLatDim, nUMLonDim, pUMTData);
					load_netcdf_2dvar(charMonthlyOBSPrcpFileName, "Prcp", iDSPMonth, nUMLatDim, nUMLonDim, pUMPData);
					//copy data from UM OBS to DS Result Array
					for(int iHighLatIndex = npMinOverlapLat; iHighLatIndex < npMaxOverlapLat; iHighLatIndex++)
					{
						for(int iHighLonIndex = npMinOverlapLon; iHighLonIndex < npMaxOverlapLon; iHighLonIndex++)
						{
							int  iCellIndex2 = iHighLatIndex*nUMLonDim + iHighLonIndex;
							float fRetriveT = *(pUMTData+iCellIndex2);
							float fRetriveP = *(pUMPData+iCellIndex2);
							*(fpDSTavgData+iCellIndex2) = fRetriveT;
							*(fpDSPrcpData+iCellIndex2) = fRetriveP;
						}
					}

				}//for(iMonth)
				printf("\n");

			}//for(lon)
		}//for(lat)
		printf("\n");
		printf("\n");
#endif

#ifdef EXPORT_BC_GRID
        //------------------------------------------
        //--------------Write BC Grid Data----------
        //------------------------------------------
		//create bc time dim
		printf("Step7: Writing Bias Correcting Dataset...");
		nGCMBCTimeDim = nGCMFUTWorkMonth;
		pGCMBCTime = (int*) malloc(sizeof(int)*nGCMBCTimeDim);
		for(int iTime=0; iTime<nGCMBCTimeDim; iTime++)
		{
			*(pGCMBCTime + iTime) = iTime + 12*(iGCMFUTWorkStartYear-1900);
		}
		//create lat and lat&lon dim
		pGCMBCLat  = pGCMFUTLat;
		pGCMBCLon  = pGCMFUTLon;
		nGCMBCLatDim  = nGCMFUTLatDim;
		nGCMBCLonDim  = nGCMFUTLonDim;
		pGCMBCData = (float*) malloc(sizeof(float) * nGCMBCLatDim*nGCMBCLonDim);

		int ncid=0, nvar_id;
		//write Tavg
		create_netcdf_3dvar(charGCMBCTavgFileName, "Tavg", nGCMBCLatDim, nGCMBCLonDim, nGCMBCTimeDim, 
			                pGCMBCData, "C", pGCMBCLat, pGCMBCLon, pGCMBCTime, 
							fGCMBCMissingValue, &ncid, &nvar_id);

		for(int iMonth = 0; iMonth < nGCMBCTimeDim; iMonth++)
		{
			pGCMBCData = grid_data_bc_Tavg[iMonth];
			save_netcdf_3dvar(pGCMBCData, nGCMBCTimeDim, nGCMBCLatDim, nGCMBCLonDim, iMonth, ncid, nvar_id);
		}
		close_netcdf_3dvar(ncid);
		//write Prcp
		create_netcdf_3dvar(charGCMBCPrcpFileName, "Prcp", nGCMBCLatDim, nGCMBCLonDim, nGCMBCTimeDim, 
			                pGCMBCData, "mm/day", pGCMBCLat, pGCMBCLon, pGCMBCTime, 
							fGCMBCMissingValue, &ncid, &nvar_id);
		for(int iMonth = 0; iMonth < nGCMBCTimeDim; iMonth++)
		{
			pGCMBCData = grid_data_bc_Prcp[iMonth];
			save_netcdf_3dvar(pGCMBCData, nGCMBCTimeDim, nGCMBCLatDim, nGCMBCLonDim, iMonth, ncid, nvar_id);
		}
		close_netcdf_3dvar(ncid);
		printf("\n");
		printf("\n");
#endif

#ifdef OUTPUT_DIAGNOSTIC
		//////////////////////////////////////////////////////////////////
		//Diagnostic Input Sequences
		//OBS, CON, FUT
		//////////////////////////////////////////////////////////////////
		sprintf(charOutputTSFileName, "%s.inputcon.csv", charInputTSFileName);
		FILE *fp_out_org;
		if( (fp_out_org = fopen(charOutputTSFileName,"wt")) == NULL)
		{
			printf("Can not open output file %s!\n", charOutputTSFileName);
			return -1;
		}
		//fprintf(fp_out_org, "TAVG_OBS,PRCP_OBS,TAVG_CON,PRCP_CON,TAVG_FUT,PRCP_FUT\n");
		
		for(int iMonth = 0; iMonth < nOBSWorkMonth; iMonth++)
		{
			for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
			{
				int iCellIndex=iRecordingCellIndex[iRecCell];
				int iLon = iDiagBCLonIndex; int iLat = iDiagBCLatIndex;
				int iCellOutIndex = iLat*nGCMFUTLonDim + iLon;
				if(iCellIndex != iCellOutIndex)
					continue;

				float *fpOBSTavgData = fOBSTSTavg[iRecCell];
				float *fpOBSPrcpData = fOBSTSPrcp[iRecCell];
				float *fpCONTavgData = fCONTSTavg[iRecCell];
				float *fpCONPrcpData = fCONTSPrcp[iRecCell];
				float fObsTAvg = *(fpOBSTavgData+iMonth);
				float fObsPrcp = *(fpOBSPrcpData+iMonth);
				float fConTAvg = *(fpCONTavgData+iMonth);
				float fConPrcp = *(fpCONPrcpData+iMonth);

				float fObsDiffTAvg = 0.0f;
				float fObsDiffPrcp = 0.0f;
				float fConDiffTAvg = 0.0f;
				float fConDiffPrcp = 0.0f;
				if(iMonth >= 12)
				{
					fObsDiffTAvg = fObsTAvg - *(fpOBSTavgData+iMonth-12);
					fObsDiffPrcp = fObsPrcp - *(fpOBSPrcpData+iMonth-12);
					fConDiffTAvg = fConTAvg - *(fpCONTavgData+iMonth-12);
					fConDiffPrcp = fConPrcp - *(fpCONPrcpData+iMonth-12);
				}
				//convert units
				fObsTAvg = fObsTAvg;
				fObsPrcp = fObsPrcp;
				fConTAvg = fConTAvg - 273.15;
				fConPrcp = fConPrcp *3600*24;
				fObsDiffTAvg = fObsDiffTAvg;
				fObsDiffPrcp = fObsDiffPrcp;
				fConDiffTAvg = fConDiffTAvg;
				fConDiffPrcp = fConDiffPrcp *3600*24;

				fprintf(fp_out_org, "%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f\n", 
									fObsTAvg, fObsPrcp, fConTAvg, fConPrcp,
									fObsDiffTAvg, fObsDiffPrcp, fConDiffTAvg, fConDiffPrcp);
			}
		}
		fclose(fp_out_org);

		sprintf(charOutputTSFileName, "%s.inputfut.csv", charInputTSFileName);
		FILE *fp_out_fut;
		if( (fp_out_fut = fopen(charOutputTSFileName,"wt")) == NULL)
		{
			printf("Can not open output file %s!\n", charOutputTSFileName);
			return -1;
		}
		//fprintf(fp_out_org, "TAVG_OBS,PRCP_OBS,TAVG_CON,PRCP_CON,TAVG_FUT,PRCP_FUT\n");
		
		for(int iMonth = 0; iMonth < nGCMFUTWorkMonth; iMonth++)
		{
			for(int iRecCell = 0; iRecCell < nRecordingCells; iRecCell++)
			{
				int iCellIndex=iRecordingCellIndex[iRecCell];
				int iLon = iDiagBCLonIndex; int iLat = iDiagBCLatIndex;
				int iCellOutIndex = iLat*nGCMFUTLonDim + iLon;
				if(iCellIndex != iCellOutIndex)
					continue;

				float *fpFUTTavgData = fFUTTSTavg[iRecCell];
				float *fpFUTPrcpData = fFUTTSPrcp[iRecCell];
				float fFutTAvg = *(fpFUTTavgData+iMonth);
				float fFutPrcp = *(fpFUTPrcpData+iMonth);

				float fFutDiffTAvg = 0.0f;
				float fFutDiffPrcp = 0.0f;
				if(iMonth >= 12)
				{
					fFutDiffTAvg = fFutTAvg - *(fpFUTTavgData+iMonth-12);
					fFutDiffPrcp = fFutPrcp - *(fpFUTPrcpData+iMonth-12);
				}
				
				//convert units
				fFutTAvg = fFutTAvg - 273.15;
				fFutPrcp = fFutPrcp *3600*24;
				fFutDiffTAvg = fFutDiffTAvg;
				fFutDiffPrcp = fFutDiffPrcp *3600*24;
				fprintf(fp_out_fut, "%5.2f,%5.2f,%5.2f,%5.2f\n", fFutTAvg, fFutPrcp, fFutDiffTAvg, fFutDiffPrcp);
			}
		}
		fclose(fp_out_fut);
#endif

#ifdef OUTPUT_DIAGNOSTIC
		//////////////////////////////////////////////////////////////////
		//Diagnostic Output BC Sequences
		//////////////////////////////////////////////////////////////////
		sprintf(charOutputTSFileName, "%s.gcmbc.csv", charInputTSFileName);

		FILE *fp_out;
		if( (fp_out = fopen(charOutputTSFileName,"wt")) == NULL)
		{
			printf("Can not open output file %s!\n", charOutputTSFileName);
			return -1;
		}
		fprintf(fp_out, "BC_TAVG,BC_PRCP\n");
		
		for(int iMonth = 12; iMonth < nGCMFUTWorkMonth; iMonth++)
		{
			float *pGCMBCTavgData = grid_data_bc_Tavg[iMonth];
			float *pGCMBCPrcpData = grid_data_bc_Prcp[iMonth];
			int iLon = iDiagBCLonIndex; int iLat = iDiagBCLatIndex;
			int iCellIndex = iLat*nGCMFUTLonDim + iLon;
			float fT = *(pGCMBCTavgData + iCellIndex);
			float fP = *(pGCMBCPrcpData + iCellIndex);
			fprintf(fp_out, "%5.2f,%5.2f\n", fT, fP);
		}
		fclose(fp_out);
#endif

#ifdef EXPORT_DS_GRID
		//------------------------------------------
        //--------------Write DS Grid Data----------
        //------------------------------------------
		//create bc time dim
		nGCMDSTimeDim = nGCMFUTWorkMonth;
		pGCMDSTime = (int*) malloc(sizeof(int)*nGCMDSTimeDim);
		for(int iTime=0; iTime<nGCMDSTimeDim; iTime++)
		{
			*(pGCMDSTime + iTime) = iTime + 12*(iGCMFUTWorkStartYear-1900);
		}
		//create lat and lat&lon dim
		pGCMDSLat     = pUMLat;
		pGCMDSLon     = pUMLon;
		nGCMDSLatDim  = nUMLatDim;
		nGCMDSLonDim  = nUMLonDim;
		pGCMDSData = (float*) malloc(sizeof(float) * nUMLatDim*nUMLonDim);
		
		//write Tavg
		//nGCMDSTimeDim = 12;
		printf("\nStep8: Writing Downscaling Dataset (Tavg)...");
		create_netcdf_3dvar(charGCMDSTavgFileName, "Tavg", nGCMDSLatDim, nGCMDSLonDim, nGCMDSTimeDim, 
			                pGCMDSData, "C", pGCMDSLat, pGCMDSLon, pGCMDSTime, 
							fGCMDSMissingValue, &ncid, &nvar_id);
		for(int iMonth = 0; iMonth < nGCMDSTimeDim; iMonth++)
		{
			printf("Writing Time Slice#=%d.\n", nGCMDSTimeDim - iMonth);
			pGCMDSData = grid_data_ds_Tavg[iMonth];
			save_netcdf_3dvar(pGCMDSData, nGCMDSTimeDim, nGCMDSLatDim, nGCMDSLonDim, iMonth, ncid, nvar_id);
		}
		close_netcdf_3dvar(ncid);
		//write Prcp
		printf("\nStep8: Writing Downscaling Dataset (Prcp)...");
		create_netcdf_3dvar(charGCMDSPrcpFileName, "Prcp", nGCMDSLatDim, nGCMDSLonDim, nGCMDSTimeDim, 
			                pGCMDSData, "mm/day", pGCMDSLat, pGCMDSLon, pGCMDSTime, 
							fGCMDSMissingValue, &ncid, &nvar_id);
		for(int iMonth = 0; iMonth < nGCMDSTimeDim; iMonth++)
		{
			printf("Writing Time Slice#=%d.\n", nGCMDSTimeDim - iMonth);
			pGCMDSData = grid_data_ds_Prcp[iMonth];
			save_netcdf_3dvar(pGCMDSData, nGCMDSTimeDim, nGCMDSLatDim, nGCMDSLonDim, iMonth, ncid, nvar_id);
		}
		close_netcdf_3dvar(ncid);
		printf("\n");
		printf("\n");
#endif

#ifdef OUTPUT_DIAGNOSTIC
		//////////////////////////////////////////////////////////////////
		//Diagnostic Output DS Sequences
		//////////////////////////////////////////////////////////////////
		sprintf(charOutputTSFileName, "%s.gcmds.csv", charInputTSFileName);

		FILE *fp_out2;
		if( (fp_out2 = fopen(charOutputTSFileName,"wt")) == NULL)
		{
			printf("Can not open output file %s!\n", charOutputTSFileName);
			return -1;
		}
		fprintf(fp_out2, "DS_TAVG,DS_PRCP\n");
		
		for(int iMonth = 12; iMonth < nGCMFUTWorkMonth; iMonth++)
		{
			float *pGCMDSTavgData = grid_data_ds_Tavg[iMonth];
			float *pGCMDSPrcpData = grid_data_ds_Prcp[iMonth];
			int iLon = iDiagDSLonIndex; int iLat = iDiagDSLatIndex;
			int iCellIndex = iLat*nUMLonDim + iLon;
			float fT = *(pGCMDSTavgData + iCellIndex);
			float fP = *(pGCMDSPrcpData + iCellIndex);
			fprintf(fp_out2, "%5.2f,%5.2f\n", fT, fP);
		}
		fclose(fp_out2);
#endif

#ifdef OUTPUT_DIAGNOSTIC
		//////////////////////////////////////////////////////////////////
		//Diagnostic Output OBS_EXPANDED Sequences
		//////////////////////////////////////////////////////////////////
		sprintf(charOutputTSFileName, "%s.obsexp.csv", charInputTSFileName);

		FILE *fp_out3;
		if( (fp_out3 = fopen(charOutputTSFileName,"wt")) == NULL)
		{
			printf("Can not open output file %s!\n", charOutputTSFileName);
			return -1;
		}
		fprintf(fp_out3, "EXP_TAVG,EXP_PRCP\n");
		
		for(int iMonth = 0; iMonth < nOBS10DTSTavg; iMonth++)
		{
			float fT = *(fOBS10DOutTavgTS + iMonth);
			float fP = *(fOBS10DOutPrcpTS + iMonth);
			fprintf(fp_out3, "%5.2f,%5.2f\n", fT, fP);
		}
		fclose(fp_out3);
#endif

        //------------------------------------------
        //--------Free Working Memory --------------
        //------------------------------------------
		for(int iMonth = 0; iMonth < num_month_bc_Tavg; iMonth++)
		{
			if(grid_data_bc_Tavg[iMonth])
				free(grid_data_bc_Tavg[iMonth]);
		}
		for(int iMonth = 0; iMonth < num_month_bc_Prcp; iMonth++)
		{
			if(grid_data_bc_Prcp[iMonth])
				free(grid_data_bc_Prcp[iMonth]);
		}
		for(int iMonth = 0; iMonth < num_month_bc_Tavg; iMonth++)
		{
			if(grid_data_dsindex_Tavg[iMonth])
				free(grid_data_dsindex_Tavg[iMonth]);
		}
		for(int iMonth = 0; iMonth < num_month_bc_Prcp; iMonth++)
		{
			if(grid_data_dsindex_Prcp[iMonth])
				free(grid_data_dsindex_Prcp[iMonth]);
		}
		for(int iMonth = 0; iMonth < num_month_ds_Tavg; iMonth++)
		{
			if(grid_data_ds_Tavg[iMonth])
				free(grid_data_ds_Tavg[iMonth]);
		}
		for(int iMonth = 0; iMonth < num_month_ds_Prcp; iMonth++)
		{
			if(grid_data_ds_Prcp[iMonth])
				free(grid_data_ds_Prcp[iMonth]);
		}
	}
	return nRetCode;
}