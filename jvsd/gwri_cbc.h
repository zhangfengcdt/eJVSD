#define MAX_DATA_LENGTH    MAX_MONTH_NUM
#include<math.h>

//function to rank the input time sequences
int get_ts_sort(float *xx, float *rr, int nx, int *iReturn)
{
	Integer exit_status=0, i, n;
	NagError fail;
	double *x=0;
	n = nx;

	INIT_FAIL(fail);
	if ( !( x = NAG_ALLOC(n, double)) )
	{
		printf("Allocation failure\n");
		*iReturn = -1;
		return -1;
	}

	for (i = 0; i < n; i++)
		x[i] = *(xx+i);


   //nag_double_sort (m01cac).
   //Quicksort of set of values of data type double
    nag_double_sort(x, (size_t) n, Nag_Ascending, &fail);

	if (fail.code != NE_NOERROR)
	{
		printf("Error from nag_ranks_and_scores (g01dhc).\n%s\n", fail.message);
		*iReturn = -1;
		return -1;
	}

	for (i = 0; i < n; i++)
		*(rr+i) = (float)x[i];

	NAG_FREE(x);
	
	*iReturn = 1;
	return 1;
}

///////////////////////////////////
//Bias Correction for a Time Series
//opt_bc=0-direct conversion;1-ECDF mapping
//opt_ecdf=0-NO 2DECDF maping;1-2DECDF joint maping
///////////////////////////////////
int ts_bc_ecdf(float* fOBSTSTavg, int nOBSTSTavg,
		   float* fOBSTSPrcp, int nOBSTSPrcp,
		   float* fCONTSTavg, int nCONTSTavg,
		   float* fCONTSPrcp, int nCONTSPrcp,
		   float* fFUTTSTavg, int nFUTTSTavg,
		   float* fFUTTSPrcp, int nFUTTSPrcp,
		   float* fBCTSTavg,  int nBCTSTavg,
		   float* fBCTSPrcp,  int nBCTSPrcp,
		   int opt_bc, int opt_ecdf)
{
	if(!fOBSTSTavg || !fOBSTSPrcp ||
	   !fCONTSTavg || !fCONTSPrcp ||
	   !fFUTTSTavg || !fFUTTSPrcp ||
	   !fBCTSTavg || !fBCTSPrcp)
	{
		printf("Fatal error: pointer to TS is null.\n");
		return -1;
	}

	if(nOBSTSTavg < 0 || nOBSTSTavg < 0 ||
	   nCONTSTavg < 0 || nCONTSPrcp < 0 ||
	   nFUTTSTavg < 0 || nFUTTSPrcp < 0 ||
	   nBCTSTavg  < 0 || nBCTSPrcp < 0)
	{
		printf("Fatal error: length of TS is negtive.\n");
		return -1;
	}

	//Input Data
	int   iDataLengthArray[DATA_COL_NUM];
	float fDataInput[DATA_COL_NUM][MAX_DATA_LENGTH];
	
	//Diff Data
	float fDataInputDiff[DATA_COL_NUM][MAX_DATA_LENGTH];
	float fDataInputDiffSort[DATA_COL_NUM][MAX_DATA_LENGTH];
	float fDiffMin[DATA_COL_NUM];
	float fDiffMax[DATA_COL_NUM];
	float fDiffMean[DATA_COL_NUM];
	float fDiffSd[DATA_COL_NUM];

	//Diff 12-Month Data
	float fDiffMonthMin[12][DATA_COL_NUM];
	float fDiffMonthMax[12][DATA_COL_NUM];
	float fDiffMonthMean[12][DATA_COL_NUM];
	float fDiffMonthSd[12][DATA_COL_NUM];

	//Output Data
	float fBCTDiffData[MAX_DATA_LENGTH];
	float fBCPDiffData[MAX_DATA_LENGTH];
	float fBCTData[MAX_DATA_LENGTH];
	float fBCPData[MAX_DATA_LENGTH];

	//////////////////////////////////////////////////////////////////
	//Step1: Reading inputs time sequences from csv files
	// The csv data file should be arranged as:
	// First line:  title
	// Second line: variables' name
	// Third line:  sequence length
	// Following lines: variables' values
	// Values orders: 
	// OBS_TAVG, OBS_PRCP, CON_TAVG, CON_PRCP, FUT_TAVG, FUT_PRCP
	//////////////////////////////////////////////////////////////////
	//colum length
	iDataLengthArray[0] = nOBSTSTavg;
	iDataLengthArray[1] = nOBSTSPrcp;
	iDataLengthArray[2] = nCONTSTavg;
	iDataLengthArray[3] = nCONTSPrcp;
	iDataLengthArray[4] = nFUTTSTavg;
	iDataLengthArray[5] = nFUTTSPrcp;

	//ts data
	for(int iLine=0; iLine<iDataLengthArray[0]; iLine++)
	{
		fDataInput[0][iLine] = *(fOBSTSTavg+iLine);
		fDataInput[1][iLine] = *(fOBSTSPrcp+iLine);
	}
	for(int iLine=0; iLine<iDataLengthArray[2]; iLine++)
	{
		fDataInput[2][iLine] = *(fCONTSTavg+iLine);
		fDataInput[3][iLine] = *(fCONTSPrcp+iLine);
	}
	for(int iLine=0; iLine<iDataLengthArray[4]; iLine++)
	{
		fDataInput[4][iLine] = *(fFUTTSTavg+iLine);
		fDataInput[5][iLine] = *(fFUTTSPrcp+iLine);
	}

	//////////////////////////////////////////////////////////////////
	//Step2: Constructing Difference Sequences
	//////////////////////////////////////////////////////////////////
	int differ_interv = 12;
	for(int i=0; i<DATA_COL_NUM; i++)
	{
		for(int iLine=0; iLine<iDataLengthArray[i]-differ_interv; iLine++)
		{
			fDataInputDiff[i][iLine] = fDataInput[i][iLine+differ_interv] 
			                          - fDataInput[i][iLine];
		}
	}

	//////////////////////////////////////////////////////////////////
	//Step3: Sorting these sequences
	//       and getting statistical summary
	//Calling NAG C Library to do this.
	//////////////////////////////////////////////////////////////////
	for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	{
		int iReturn = 1;
		get_ts_sort(&fDataInputDiff[iSeq][0], &fDataInputDiffSort[iSeq][0], iDataLengthArray[iSeq] - differ_interv, &iReturn);
		if(iReturn < 0)
		{
			printf("Can not sort diff data!\n");
			return -1;
		}
	}

	for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	{
		double xsd, xskew, xkurt, wsum, xmean, xmax, xmin;
		double x[MAX_DATA_LENGTH];
		Integer n, nvalid;
		static NagError fail;

		n = iDataLengthArray[iSeq] - differ_interv;
		for (int i = 0; i < iDataLengthArray[iSeq] - differ_interv; i++)
			x[i] = (fDataInputDiff[iSeq][i]);


		nag_summary_stats_1var(n, x, (double *)0, &nvalid, &xmean, &xsd, &xskew, &xkurt,
			&xmin, &xmax, &wsum, &fail);

		if (fail.code != NE_NOERROR)
		{
			printf("Error from nag_ranks_and_scores (g01dhc).\n%s\n", fail.message);
			return -1;
		}

		fDiffMin[iSeq]  = (float)xmin;
		fDiffMax[iSeq]  = (float)xmax;
		fDiffMean[iSeq] = (float)xmean;
		fDiffSd[iSeq]   = (float)xsd;
	}

	//for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	//{
	//	double xsd, xskew, xkurt, wsum, xmean, xmax, xmin;
	//	double x[MAX_DATA_LENGTH];
	//	Integer n, nvalid;
	//	static NagError fail;

	//	for(int iMonth=0; iMonth<12; iMonth++)
	//	{
	//		n = (iDataLengthArray[iSeq] - differ_interv)/12;
	//		for (int i = 0; i < (iDataLengthArray[iSeq] - differ_interv)/12; i++)
	//			x[i] = (fDataInputDiff[iSeq][i*12+iMonth]);

	//		nag_summary_stats_1var(n, x, (double *)0, &nvalid, &xmean, &xsd, &xskew, &xkurt,
	//			&xmin, &xmax, &wsum, &fail);

	//		if (fail.code != NE_NOERROR)
	//		{
	//			printf("Error from nag_ranks_and_scores (g01dhc).\n%s\n", fail.message);
	//			return -1;
	//		}

	//		fDiffMonthMin[iMonth][iSeq]  = (float)xmin;
	//		fDiffMonthMax[iMonth][iSeq]  = (float)xmax;
	//		fDiffMonthMean[iMonth][iSeq] = (float)xmean;
	//		fDiffMonthSd[iMonth][iSeq]   = (float)xsd;
	//	}
	//}

	//////////////////////////////////////////////////////////////////
	//Step4: Bias Correction for future sequences
	//       by using marginal cdf and joint cdf
	//////////////////////////////////////////////////////////////////
	for(int iDiffIndex=0; iDiffIndex < iDataLengthArray[4] - differ_interv; iDiffIndex++)
	{
		//retrive data
		float fTavgDiff = fDataInputDiff[4][iDiffIndex];
		float fPrcpDiff = fDataInputDiff[5][iDiffIndex];

		float std_TavgDiff_OBS = fDiffSd[0];
		float std_PrcpDiff_OBS = fDiffSd[1];
		float std_TavgDiff_CON = fDiffSd[2];
		float std_PrcpDiff_CON = fDiffSd[3];
		float std_TavgDiff_FUT = fDiffSd[4];
		float std_PrcpDiff_FUT = fDiffSd[5];

		if(opt_bc == 1)
		{
			//Marginal ECDF Mapping Method
			//find marginal frequency of TavgDiff on Control Run
			float fTavgMarginalFreq, fPrcpMarginalFreq;
			int iData;
			fTavgMarginalFreq = 1.0f;
			for(iData=0; iData < iDataLengthArray[2] - differ_interv; iData++)
			{
				if(fDataInputDiffSort[2][iData] >= fTavgDiff)
				{
					if(iData == 0) //out of minimum limit
					{
						fTavgMarginalFreq = 0.0f;
						break;
					}
					else
					{
						float fRank = iData - 1 + (fTavgDiff - fDataInputDiffSort[2][iData-1])/
												  (fDataInputDiffSort[2][iData]-fDataInputDiffSort[2][iData-1]);
						fTavgMarginalFreq = fRank/(iDataLengthArray[2] - differ_interv);
						break;
					}
				}
			}
			if(iData == iDataLengthArray[2] - differ_interv -1) //out of maximum limit
			{
				fTavgMarginalFreq = 1.0f;
			}
			//find marginal frequency of PrcpDiff on Control Run
			fPrcpMarginalFreq = 1.0f;
			for(iData=0; iData < iDataLengthArray[3] - differ_interv; iData++)
			{
				if(fDataInputDiffSort[3][iData] > fPrcpDiff)
				{
					if(iData == 0) //out of minimum limit
					{
						fPrcpMarginalFreq = 0.0f;
						break;
					}
					else
					{
						float fRank = iData - 1 + (fPrcpDiff - fDataInputDiffSort[3][iData-1])/
												  (fDataInputDiffSort[3][iData]-fDataInputDiffSort[3][iData-1]);
						fPrcpMarginalFreq = fRank/(iDataLengthArray[3] - differ_interv);
						break;
					}
				}
			}
			if(iData == iDataLengthArray[2] - differ_interv) //out of maximum limit
			{
				fPrcpMarginalFreq = 1.0f;
			}

			//interpolate TavgDiff on OBS using marginal cdf
			if(fTavgMarginalFreq == 0.0f)
				fBCTDiffData[iDiffIndex] = fDataInputDiffSort[0][0];
			else if(fTavgMarginalFreq == 1.0f)
				fBCTDiffData[iDiffIndex] = fDataInputDiffSort[0][iDataLengthArray[0] - differ_interv - 1];
			else
			{
				int iDataLength = iDataLengthArray[0] - differ_interv;
				int iLeftIndex  = int(fTavgMarginalFreq*iDataLength);

				if(iLeftIndex >=  iDataLength -1)
					fBCTDiffData[iDiffIndex] = fDataInputDiffSort[0][iDataLength - 1];
				else
				{
					float fDelt     = fDataInputDiffSort[0][iLeftIndex+1] - fDataInputDiffSort[0][iLeftIndex];
					float fDeltFreq = (float)(iLeftIndex+1)/(float)iDataLength;
					float fDiffFreq = (float)iLeftIndex/(float)iDataLength;
					fBCTDiffData[iDiffIndex] = fDataInputDiffSort[0][iLeftIndex]
									+ (fTavgMarginalFreq -  fDiffFreq)/(fDeltFreq - fDiffFreq)* fDelt;
				}
			}

			//interpolate PrcpDiff on OBS using marginal cdf
			if(fPrcpMarginalFreq == 0.0f)
				fBCPDiffData[iDiffIndex] = fDataInputDiffSort[1][0];
			else if(fPrcpMarginalFreq == 1.0f)
				fBCPDiffData[iDiffIndex] = fDataInputDiffSort[1][iDataLengthArray[1] - differ_interv - 1];
			else
			{
				int iDataLength = iDataLengthArray[1] - differ_interv;
				int iLeftIndex  = int(fPrcpMarginalFreq*iDataLength);

				if(iLeftIndex >=  iDataLength -1)
					fBCPDiffData[iDiffIndex] = fDataInputDiffSort[1][iDataLength - 1];
				else
				{
					float fDelt     = fDataInputDiffSort[1][iLeftIndex+1] - fDataInputDiffSort[1][iLeftIndex];
					float fDeltFreq = (float)(iLeftIndex+1)/(float)iDataLength;
					float fDiffFreq = (float)iLeftIndex/(float)iDataLength;
					fBCPDiffData[iDiffIndex] = fDataInputDiffSort[1][iLeftIndex]
									+ (fPrcpMarginalFreq -  fDiffFreq)/(fDeltFreq - fDiffFreq)* fDelt;
				}
			}

			if(opt_ecdf == 1)
			{
				float fFreqDiffCon, fFreqDiffObs;

				//find the 2D cdf value on control run
				int iCount=0;
				for(int i = 0; i < iDataLengthArray[2] - differ_interv; i++)
				{
					float fTDiff = fDataInputDiff[2][i];
					float fPDiff = fDataInputDiff[3][i];

					if(fTDiff < fTavgDiff && fPDiff < fPrcpDiff)
						iCount++;
				}
				fFreqDiffCon = (float)iCount / (float) (iDataLengthArray[2] - differ_interv);

				//find the 2D cdf value on observation
				iCount=0;
				for(int i = 0; i < iDataLengthArray[0] - differ_interv; i++)
				{
					float fTDiff = fDataInputDiff[0][i];
					float fPDiff = fDataInputDiff[1][i];

					if(fTDiff < fBCTDiffData[iDiffIndex] && fPDiff < fBCPDiffData[iDiffIndex])
						iCount++;
				}
				fFreqDiffObs = (float)iCount / (float) (iDataLengthArray[0] - differ_interv);

				float fFreqDiff = fFreqDiffCon - fFreqDiffObs;
				if(fFreqDiff > -LMT_JTMG_CDFDIFF && fFreqDiff < LMT_JTMG_CDFDIFF)
					continue;

				if(fFreqDiffCon < fFreqDiffObs)
				{
					for(int iBin=0; iBin<=NUM_BIN; iBin++)
					{
						float fTDiffLimit = fBCTDiffData[iDiffIndex] - (float)iBin/(float)NUM_BIN * (fBCTDiffData[iDiffIndex] - fDiffMin[0]);
						float fPDiffLimit = fBCPDiffData[iDiffIndex] - (float)iBin/(float)NUM_BIN * (fBCPDiffData[iDiffIndex] - fDiffMin[1]);

						int iCount=0;
						for(int i = 0; i < iDataLengthArray[0] - differ_interv; i++)
						{
							float fTDiff = fDataInputDiff[0][i];
							float fPDiff = fDataInputDiff[1][i];

							if(fTDiff < fTDiffLimit && fPDiff < fPDiffLimit)
								iCount++;
						}
						float fFreqDiffObs2 = (float)iCount / (float) (iDataLengthArray[0] - differ_interv);
						float fFreqDiffInter = fFreqDiffCon - fFreqDiffObs2;
						if(fFreqDiffInter > -LMT_JTMG_CDFDIFF && fFreqDiffInter < LMT_JTMG_CDFDIFF)
						{
							fBCTDiffData[iDiffIndex] = fTDiffLimit;
							fBCPDiffData[iDiffIndex] = fPDiffLimit;
							break;
						}
					}
				}
				else
				{
					for(int iBin=0; iBin<=NUM_BIN; iBin++)
					{
						float fTDiffLimit = fBCTDiffData[iDiffIndex] + (float)iBin/(float)NUM_BIN * (fDiffMax[0] - fBCTDiffData[iDiffIndex]);
						float fPDiffLimit = fBCPDiffData[iDiffIndex] + (float)iBin/(float)NUM_BIN * (fDiffMax[1] - fBCPDiffData[iDiffIndex]);

						int iCount=0;
						for(int i = 0; i < iDataLengthArray[0] - differ_interv; i++)
						{
							float fTDiff = fDataInputDiff[0][i];
							float fPDiff = fDataInputDiff[1][i];

							if(fTDiff < fTDiffLimit && fPDiff < fPDiffLimit)
								iCount++;
						}
						float fFreqDiffObs2 = (float)iCount / (float) (iDataLengthArray[0] - differ_interv);
						float fFreqDiffInter = fFreqDiffCon - fFreqDiffObs2;
						if(fFreqDiffInter > -LMT_JTMG_CDFDIFF && fFreqDiffInter < LMT_JTMG_CDFDIFF)
						{
							fBCTDiffData[iDiffIndex] = fTDiffLimit;
							fBCPDiffData[iDiffIndex] = fPDiffLimit;
							break;
						}
					}
				}
			}//if(opt_ecdf == 1)

		}
		else //(opt_bc != 1)
		{
			float std_TavgDiff_OBS = fDiffMonthSd[iDiffIndex%12][0];
			float std_PrcpDiff_OBS = fDiffMonthSd[iDiffIndex%12][1];
			float std_TavgDiff_CON = fDiffMonthSd[iDiffIndex%12][2];
			float std_PrcpDiff_CON = fDiffMonthSd[iDiffIndex%12][3];
			float std_TavgDiff_FUT = fDiffMonthSd[iDiffIndex%12][4];
			float std_PrcpDiff_FUT = fDiffMonthSd[iDiffIndex%12][5];

			float a1  = std_TavgDiff_OBS/std_TavgDiff_CON;
			float b1  = std_PrcpDiff_OBS/std_PrcpDiff_CON;
			float ro1 = -0.1;
			float ro2 = -0.1;

			float std_TavgDiff_BC = a1*std_TavgDiff_FUT;
			float std_PrcpDiff_BC = b1*std_PrcpDiff_FUT;

			if(opt_ecdf == 1)
			{
				//equations need to be revised
				fBCTDiffData[iDiffIndex] = std_TavgDiff_BC/std_TavgDiff_FUT * fTavgDiff 
											+ std_TavgDiff_FUT/std_PrcpDiff_FUT*ro1*fPrcpDiff;
				fBCPDiffData[iDiffIndex] = std_PrcpDiff_BC/std_PrcpDiff_FUT * fPrcpDiff
											+ std_PrcpDiff_FUT/std_TavgDiff_FUT*ro2*fTavgDiff;
			}
			else
			{
				fBCTDiffData[iDiffIndex] = fTavgDiff/std_TavgDiff_FUT * std_TavgDiff_BC;
				fBCPDiffData[iDiffIndex] = fPrcpDiff/std_PrcpDiff_FUT * std_PrcpDiff_BC;
			}
		}


	}//for(iDiffIndex)

	//////////////////////////////////////////////////////////////////
	//Step5: Reconsturcting Tavg and Prcp from
	//       the bias corrected differces sequences
	//////////////////////////////////////////////////////////////////
	//calculate climatology of obs/con/fut
	float fClimat[DATA_COL_NUM][12];
	for(int iCol=0; iCol<DATA_COL_NUM; iCol++)
	{
		for(int iMonth=0; iMonth<12; iMonth++)
		{
			float fTotalSum = 0.0f;
			for(int iData=0+iMonth; iData<iDataLengthArray[iCol]; iData+=12)
			{
				fTotalSum += fDataInput[iCol][iData];
			}
			fClimat[iCol][iMonth] = fTotalSum / (float) (iDataLengthArray[iCol]/12);
		}
	}

	//pure incremental series
	for(int iData=0; iData<iDataLengthArray[4]; iData++)
	{
		if(iData<12)
		{
			fBCTData[iData] = 0.0f;
			fBCPData[iData] = 0.0f;
		}
		else
		{
			fBCTData[iData] = fBCTData[iData-12] + fBCTDiffData[iData-12];
			fBCPData[iData] = fBCPData[iData-12] + fBCPDiffData[iData-12];
			//fBCTData[iData] = fBCTData[iData-12] + fDataInputDiff[4][iData-12];
			//fBCPData[iData] = fBCPData[iData-12] + fDataInputDiff[5][iData-12];
		}
	}

	//calculate climatology of pure incremental series
	float fClimatIncTavg[12], fClimatIncPrcp[12];
	for(int iMonth=0; iMonth<12; iMonth++)
	{
		float fTotalSumTavg = 0.0f;
		float fTotalSumPrcp = 0.0f;
		for(int iData=0+iMonth+12; iData<iDataLengthArray[4]; iData+=12)
		{
			fTotalSumTavg += fBCTData[iData];
			fTotalSumPrcp += fBCPData[iData];
		}
		fClimatIncTavg[iMonth] = fTotalSumTavg / (float) (iDataLengthArray[4]/12-1);
		fClimatIncPrcp[iMonth] = fTotalSumPrcp / (float) (iDataLengthArray[4]/12-1);
	}

	//adding bias into pure incremental series
	for(int iData=0; iData<iDataLengthArray[4]; iData++)
	{
		if(iData<12)
		{
			fBCTData[iData] = fClimat[4][iData] - (fClimat[2][iData] - fClimat[0][iData]) - fClimatIncTavg[iData];
			fBCPData[iData] = fClimat[5][iData] - (fClimat[3][iData] - fClimat[1][iData]) - fClimatIncPrcp[iData];
		}
		else
		{
			fBCTData[iData] = fBCTData[iData] + fBCTData[iData%12];
			fBCPData[iData] = fBCPData[iData] + fBCPData[iData%12];
		}
	}

	//make sure the precipitation is positive
	for(int i=0; i < iDataLengthArray[4]; i++)
	{
		if(fBCPData[i] < 0) 
		{
			//negtive extreme value processing
			fBCPData[i] = 1/(1-fBCPData[i]);
		}
	}

	////////////////////////////////////////////////////////////////////
	////Step6: Writting BC Sequences
	////////////////////////////////////////////////////////////////////
	for(int i=0; i < iDataLengthArray[4]; i++)
	{
		* (fBCTSTavg+i) = fBCTData[i];
		* (fBCTSPrcp+i) = fBCPData[i];
	}

	return 0;
}


///////////////////////////////////
//Bias Correction for a Time Series
//Using Anomany Sequences
//opt_bc=0-direct conversion;1-ECDF mapping
//opt_ecdf=0-NO 2DECDF maping;1-2DECDF joint maping
///////////////////////////////////
int ts_bc_ecdf2(float* fOBSTSTavg, int nOBSTSTavg,
		   float* fOBSTSPrcp, int nOBSTSPrcp,
		   float* fCONTSTavg, int nCONTSTavg,
		   float* fCONTSPrcp, int nCONTSPrcp,
		   float* fFUTTSTavg, int nFUTTSTavg,
		   float* fFUTTSPrcp, int nFUTTSPrcp,
		   float* fBCTSTavg,  int nBCTSTavg,
		   float* fBCTSPrcp,  int nBCTSPrcp,
		   int opt_bc, int opt_ecdf)
{
	if(!fOBSTSTavg || !fOBSTSPrcp ||
	   !fCONTSTavg || !fCONTSPrcp ||
	   !fFUTTSTavg || !fFUTTSPrcp ||
	   !fBCTSTavg || !fBCTSPrcp)
	{
		printf("Fatal error: pointer to TS is null.\n");
		return -1;
	}

	if(nOBSTSTavg < 0 || nOBSTSTavg < 0 ||
	   nCONTSTavg < 0 || nCONTSPrcp < 0 ||
	   nFUTTSTavg < 0 || nFUTTSPrcp < 0 ||
	   nBCTSTavg  < 0 || nBCTSPrcp < 0)
	{
		printf("Fatal error: length of TS is negtive.\n");
		return -1;
	}

	//Input Data
	int   iDataLengthArray[DATA_COL_NUM];
	float fDataInput[DATA_COL_NUM][MAX_DATA_LENGTH];
	
	//Diff Data
	float fDataInputDiff[DATA_COL_NUM][MAX_DATA_LENGTH];
	float fDataInputDiffSort[DATA_COL_NUM][MAX_DATA_LENGTH];
	float fDataInputDiffMonthly[12][DATA_COL_NUM][MAX_DATA_LENGTH/12];
	float fDataInputDiffMonthlySort[12][DATA_COL_NUM][MAX_DATA_LENGTH/12];
	float fDiffMin[DATA_COL_NUM];
	float fDiffMax[DATA_COL_NUM];
	float fDiffMean[DATA_COL_NUM];
	float fDiffSd[DATA_COL_NUM];

	//Diff 12-Month Data
	float fDiffMonthMin[12][DATA_COL_NUM];
	float fDiffMonthMax[12][DATA_COL_NUM];
	float fDiffMonthMean[12][DATA_COL_NUM];
	float fDiffMonthSd[12][DATA_COL_NUM];

	//Output Data
	float fBCTDiffData[MAX_DATA_LENGTH];
	float fBCPDiffData[MAX_DATA_LENGTH];
	float fBCTData[MAX_DATA_LENGTH];
	float fBCPData[MAX_DATA_LENGTH];

	//////////////////////////////////////////////////////////////////
	//Step1: Reading inputs time sequences from csv files
	// The csv data file should be arranged as:
	// First line:  title
	// Second line: variables' name
	// Third line:  sequence length
	// Following lines: variables' values
	// Values orders: 
	// OBS_TAVG, OBS_PRCP, CON_TAVG, CON_PRCP, FUT_TAVG, FUT_PRCP
	//////////////////////////////////////////////////////////////////
	//colum length
	iDataLengthArray[0] = nOBSTSTavg;
	iDataLengthArray[1] = nOBSTSPrcp;
	iDataLengthArray[2] = nCONTSTavg;
	iDataLengthArray[3] = nCONTSPrcp;
	iDataLengthArray[4] = nFUTTSTavg;
	iDataLengthArray[5] = nFUTTSPrcp;

	//ts data
	for(int iLine=0; iLine<iDataLengthArray[0]; iLine++)
	{
		fDataInput[0][iLine] = *(fOBSTSTavg+iLine);
		fDataInput[1][iLine] = *(fOBSTSPrcp+iLine);
	}
	for(int iLine=0; iLine<iDataLengthArray[2]; iLine++)
	{
		fDataInput[2][iLine] = *(fCONTSTavg+iLine);
		fDataInput[3][iLine] = *(fCONTSPrcp+iLine);
	}
	for(int iLine=0; iLine<iDataLengthArray[4]; iLine++)
	{
		fDataInput[4][iLine] = *(fFUTTSTavg+iLine);
		fDataInput[5][iLine] = *(fFUTTSPrcp+iLine);
	}

	//////////////////////////////////////////////////////////////////
	//Step2: Constructing Difference Sequences
	//////////////////////////////////////////////////////////////////
	int differ_interv = 0;
	float temp_mean[DATA_COL_NUM][12];
	for(int iCol=0; iCol<DATA_COL_NUM; iCol++)
	{
		for(int iMonth=0; iMonth<12; iMonth++)
		{
			float fTotalSum = 0.0f;
			for(int iData=0+iMonth; iData<iDataLengthArray[iCol]; iData+=12)
			{
				fTotalSum += fDataInput[iCol][iData];
			}
			temp_mean[iCol][iMonth] = fTotalSum / (float) (iDataLengthArray[iCol]/12);
		}
	}

	for(int iLine=0; iLine<iDataLengthArray[0]; iLine++)
	{
			fDataInputDiff[0][iLine] = fDataInput[0][iLine] - temp_mean[0][iLine%12];
			fDataInputDiff[2][iLine] = fDataInput[2][iLine] - temp_mean[2][iLine%12];
			fDataInputDiff[4][iLine] = fDataInput[4][iLine] - temp_mean[4][iLine%12];
			//fDataInputDiff[2][iLine] = fDataInput[0][iLine] - temp_mean[0][iLine%12];
			//fDataInputDiff[4][iLine] = fDataInput[0][iLine] - temp_mean[0][iLine%12];
	}
	for(int iLine=0; iLine<iDataLengthArray[0]-differ_interv; iLine++)
	{
			fDataInputDiff[1][iLine] = fDataInput[1][iLine] - temp_mean[1][iLine%12];
			fDataInputDiff[3][iLine] = fDataInput[3][iLine] - temp_mean[3][iLine%12];
			fDataInputDiff[5][iLine] = fDataInput[5][iLine] - temp_mean[5][iLine%12];
			//fDataInputDiff[3][iLine] = fDataInput[1][iLine] - temp_mean[1][iLine%12];
			//fDataInputDiff[5][iLine] = fDataInput[1][iLine] - temp_mean[1][iLine%12];
	}

	//////////////////////////////////////////////////////////////////
	//Step2.2: Generating Monthly Diff Array
	//////////////////////////////////////////////////////////////////
	for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	{
		for(int iRec=0; iRec<iDataLengthArray[iSeq] - differ_interv; iRec++)
		{
			int iMonth = iRec%12;
			int iYear  = (int)(iRec/12);
			fDataInputDiffMonthly[iMonth][iSeq][iYear] = fDataInputDiff[iSeq][iRec];
		}
	}


	//////////////////////////////////////////////////////////////////
	//Step3: Sorting these sequences
	//       and getting statistical summary
	//Calling NAG C Library to do this.
	//////////////////////////////////////////////////////////////////
	for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	{
		int iReturn = 1;
		get_ts_sort(&fDataInputDiff[iSeq][0], &fDataInputDiffSort[iSeq][0], iDataLengthArray[iSeq] - differ_interv, &iReturn);
		if(iReturn < 0)
		{
			printf("Can not sort diff data!\n");
			return -1;
		}
		for(int iMonth=0; iMonth<12; iMonth++)
		{
			int iReturn = 1;
			get_ts_sort(&fDataInputDiffMonthly[iMonth][iSeq][0], &fDataInputDiffMonthlySort[iMonth][iSeq][0], (iDataLengthArray[iSeq] - differ_interv)/12, &iReturn);
			if(iReturn < 0)
			{
				printf("Can not sort diff data!\n");
				return -1;
			}
		}
	}

	for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	{
		double xsd, xskew, xkurt, wsum, xmean, xmax, xmin;
		double x[MAX_DATA_LENGTH];
		Integer n, nvalid;
		static NagError fail;

		n = iDataLengthArray[iSeq] - differ_interv;
		for (int i = 0; i < iDataLengthArray[iSeq] - differ_interv; i++)
			x[i] = (fDataInputDiff[iSeq][i]);


		nag_summary_stats_1var(n, x, (double *)0, &nvalid, &xmean, &xsd, &xskew, &xkurt,
			&xmin, &xmax, &wsum, &fail);

		if (fail.code != NE_NOERROR)
		{
			printf("Error from nag_ranks_and_scores (g01dhc).\n%s\n", fail.message);
			return -1;
		}

		fDiffMin[iSeq]  = (float)xmin;
		fDiffMax[iSeq]  = (float)xmax;
		fDiffMean[iSeq] = (float)xmean;
		fDiffSd[iSeq]   = (float)xsd;
	}

	for(int iSeq=0; iSeq<DATA_COL_NUM; iSeq++)
	{
		double xsd, xskew, xkurt, wsum, xmean, xmax, xmin;
		double x[MAX_DATA_LENGTH];
		Integer n, nvalid;
		static NagError fail;

		for(int iMonth=0; iMonth<12; iMonth++)
		{
			n = (iDataLengthArray[iSeq] - differ_interv)/12;
			for (int i = 0; i < (iDataLengthArray[iSeq] - differ_interv)/12; i++)
				x[i] = (fDataInputDiff[iSeq][i*12+iMonth]);

			nag_summary_stats_1var(n, x, (double *)0, &nvalid, &xmean, &xsd, &xskew, &xkurt,
				&xmin, &xmax, &wsum, &fail);

			if (fail.code != NE_NOERROR)
			{
				printf("Error from nag_ranks_and_scores (g01dhc).\n%s\n", fail.message);
				return -1;
			}

			fDiffMonthMin[iMonth][iSeq]  = (float)xmin;
			fDiffMonthMax[iMonth][iSeq]  = (float)xmax;
			fDiffMonthMean[iMonth][iSeq] = (float)xmean;
			fDiffMonthSd[iMonth][iSeq]   = (float)xsd;
		}
	}

	//////////////////////////////////////////////////////////////////
	//Step4: Bias Correction for future sequences
	//       by using marginal cdf and joint cdf
	//////////////////////////////////////////////////////////////////
	for(int iDiffIndex=0; iDiffIndex < iDataLengthArray[4] - differ_interv; iDiffIndex++)
	{
		//retrive data
		float fTavgDiff = fDataInputDiff[4][iDiffIndex];
		float fPrcpDiff = fDataInputDiff[5][iDiffIndex];

		int iMonth = iDiffIndex%12;

		float std_TavgDiff_OBS = fDiffSd[0];
		float std_PrcpDiff_OBS = fDiffSd[1];
		float std_TavgDiff_CON = fDiffSd[2];
		float std_PrcpDiff_CON = fDiffSd[3];
		float std_TavgDiff_FUT = fDiffSd[4];
		float std_PrcpDiff_FUT = fDiffSd[5];

		if(opt_bc == 1)
		{
			
			//Marginal ECDF Mapping Method
			//find marginal frequency of TavgDiff on Control Run
			float fTavgMarginalFreq, fPrcpMarginalFreq;
			int iData;
			fTavgMarginalFreq = 1.0f;
			for(iData=0; iData < (iDataLengthArray[2] - differ_interv)/12; iData++)
			{
				if(fDataInputDiffMonthlySort[iMonth][2][iData] >= fTavgDiff)
				{
					if(iData == 0) //out of minimum limit
					{
						fTavgMarginalFreq = 0.0f;
						break;
					}
					else
					{
						float fRank = iData - 1 + (fTavgDiff - fDataInputDiffMonthlySort[iMonth][2][iData-1])/
												  (fDataInputDiffMonthlySort[iMonth][2][iData]-fDataInputDiffMonthlySort[iMonth][2][iData-1]);
						fTavgMarginalFreq = fRank/((iDataLengthArray[2] - differ_interv)/12);
						break;
					}
				}
			}
			if(iData == (iDataLengthArray[2] - differ_interv)/12 -1) //out of maximum limit
			{
				fTavgMarginalFreq = 1.0f;
			}
			//find marginal frequency of PrcpDiff on Control Run
			fPrcpMarginalFreq = 1.0f;
			for(iData=0; iData < (iDataLengthArray[3] - differ_interv)/12; iData++)
			{
				if(fDataInputDiffMonthlySort[iMonth][3][iData] > fPrcpDiff)
				{
					if(iData == 0) //out of minimum limit
					{
						fPrcpMarginalFreq = 0.0f;
						break;
					}
					else
					{
						float fRank = iData - 1 + (fPrcpDiff - fDataInputDiffMonthlySort[iMonth][3][iData-1])/
												  (fDataInputDiffMonthlySort[iMonth][3][iData]-fDataInputDiffMonthlySort[iMonth][3][iData-1]);
						fPrcpMarginalFreq = fRank/((iDataLengthArray[3] - differ_interv)/12);
						break;
					}
				}
			}
			if(iData == (iDataLengthArray[2] - differ_interv)/2 -1) //out of maximum limit
			{
				fPrcpMarginalFreq = 1.0f;
			}

			//interpolate TavgDiff on OBS using marginal cdf
			if(fTavgMarginalFreq == 0.0f)
				fBCTDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][0][0];
			else if(fTavgMarginalFreq == 1.0f)
				fBCTDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][0][(iDataLengthArray[0] - differ_interv)/12 - 1];
			else
			{
				int iDataLength = (iDataLengthArray[0] - differ_interv)/12;
				int iLeftIndex  = int(fTavgMarginalFreq*iDataLength);

				if(iLeftIndex >=  iDataLength -1)
					fBCTDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][0][iDataLength - 1];
				else
				{
					float fDelt     = fDataInputDiffMonthlySort[iMonth][0][iLeftIndex+1] - fDataInputDiffMonthlySort[iMonth][0][iLeftIndex];
					float fDeltFreq = (float)(iLeftIndex+1)/(float)iDataLength;
					float fDiffFreq = (float)iLeftIndex/(float)iDataLength;
					fBCTDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][0][iLeftIndex]
									+ (fTavgMarginalFreq -  fDiffFreq)/(fDeltFreq - fDiffFreq)* fDelt;
				}
			}

			//interpolate PrcpDiff on OBS using marginal cdf
			if(fPrcpMarginalFreq == 0.0f)
				fBCPDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][1][0];
			else if(fPrcpMarginalFreq == 1.0f)
				fBCPDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][1][(iDataLengthArray[1] - differ_interv)/12 - 1];
			else
			{
				int iDataLength = (iDataLengthArray[1] - differ_interv)/12;
				int iLeftIndex  = int(fPrcpMarginalFreq*iDataLength);

				if(iLeftIndex >=  iDataLength -1)
					fBCPDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][1][iDataLength - 1];
				else
				{
					float fDelt     = fDataInputDiffMonthlySort[iMonth][1][iLeftIndex+1] - fDataInputDiffMonthlySort[iMonth][1][iLeftIndex];
					float fDeltFreq = (float)(iLeftIndex+1)/(float)iDataLength;
					float fDiffFreq = (float)iLeftIndex/(float)iDataLength;
					fBCPDiffData[iDiffIndex] = fDataInputDiffMonthlySort[iMonth][1][iLeftIndex]
									+ (fPrcpMarginalFreq -  fDiffFreq)/(fDeltFreq - fDiffFreq)* fDelt;
				}
			}

			if(opt_ecdf == 1)
			{
				float fFreqDiffCon, fFreqDiffObs;

				//find the 2D cdf value on control run
				int iCount=0;
				for(int i = 0; i < iDataLengthArray[2] - differ_interv; i++)
				{
					float fTDiff = fDataInputDiff[2][i];
					float fPDiff = fDataInputDiff[3][i];

					if(fTDiff < fTavgDiff && fPDiff < fPrcpDiff)
						iCount++;
				}
				fFreqDiffCon = (float)iCount / (float) (iDataLengthArray[2] - differ_interv);

				//find the 2D cdf value on observation
				iCount=0;
				for(int i = 0; i < iDataLengthArray[0] - differ_interv; i++)
				{
					float fTDiff = fDataInputDiff[0][i];
					float fPDiff = fDataInputDiff[1][i];

					if(fTDiff < fBCTDiffData[iDiffIndex] && fPDiff < fBCPDiffData[iDiffIndex])
						iCount++;
				}
				fFreqDiffObs = (float)iCount / (float) (iDataLengthArray[0] - differ_interv);

				float fFreqDiff = fFreqDiffCon - fFreqDiffObs;
				if(fFreqDiff > -LMT_JTMG_CDFDIFF && fFreqDiff < LMT_JTMG_CDFDIFF)
					continue;

				if(fFreqDiffCon < fFreqDiffObs)
				{
					for(int iBin=0; iBin<=NUM_BIN; iBin++)
					{
						float fTDiffLimit = fBCTDiffData[iDiffIndex] - (float)iBin/(float)NUM_BIN * (fBCTDiffData[iDiffIndex] - fDiffMin[0]);
						float fPDiffLimit = fBCPDiffData[iDiffIndex] - (float)iBin/(float)NUM_BIN * (fBCPDiffData[iDiffIndex] - fDiffMin[1]);

						int iCount=0;
						for(int i = 0; i < iDataLengthArray[0] - differ_interv; i++)
						{
							float fTDiff = fDataInputDiff[0][i];
							float fPDiff = fDataInputDiff[1][i];

							if(fTDiff < fTDiffLimit && fPDiff < fPDiffLimit)
								iCount++;
						}
						float fFreqDiffObs2 = (float)iCount / (float) (iDataLengthArray[0] - differ_interv);
						float fFreqDiffInter = fFreqDiffCon - fFreqDiffObs2;
						if(fFreqDiffInter > -LMT_JTMG_CDFDIFF && fFreqDiffInter < LMT_JTMG_CDFDIFF)
						{
							fBCTDiffData[iDiffIndex] = fTDiffLimit;
							fBCPDiffData[iDiffIndex] = fPDiffLimit;
							break;
						}
					}
				}
				else
				{
					for(int iBin=0; iBin<=NUM_BIN; iBin++)
					{
						float fTDiffLimit = fBCTDiffData[iDiffIndex] + (float)iBin/(float)NUM_BIN * (fDiffMax[0] - fBCTDiffData[iDiffIndex]);
						float fPDiffLimit = fBCPDiffData[iDiffIndex] + (float)iBin/(float)NUM_BIN * (fDiffMax[1] - fBCPDiffData[iDiffIndex]);

						int iCount=0;
						for(int i = 0; i < iDataLengthArray[0] - differ_interv; i++)
						{
							float fTDiff = fDataInputDiff[0][i];
							float fPDiff = fDataInputDiff[1][i];

							if(fTDiff < fTDiffLimit && fPDiff < fPDiffLimit)
								iCount++;
						}
						float fFreqDiffObs2 = (float)iCount / (float) (iDataLengthArray[0] - differ_interv);
						float fFreqDiffInter = fFreqDiffCon - fFreqDiffObs2;
						if(fFreqDiffInter > -LMT_JTMG_CDFDIFF && fFreqDiffInter < LMT_JTMG_CDFDIFF)
						{
							fBCTDiffData[iDiffIndex] = fTDiffLimit;
							fBCPDiffData[iDiffIndex] = fPDiffLimit;
							break;
						}
					}
				}
			}//if(opt_ecdf == 1)

		}
		else //(opt_bc != 1)
		{
			float std_TavgDiff_OBS = fDiffMonthSd[iDiffIndex%12][0];
			float std_PrcpDiff_OBS = fDiffMonthSd[iDiffIndex%12][1];
			float std_TavgDiff_CON = fDiffMonthSd[iDiffIndex%12][2];
			float std_PrcpDiff_CON = fDiffMonthSd[iDiffIndex%12][3];
			float std_TavgDiff_FUT = fDiffMonthSd[iDiffIndex%12][4];
			float std_PrcpDiff_FUT = fDiffMonthSd[iDiffIndex%12][5];

			float a1  = std_TavgDiff_OBS/std_TavgDiff_CON;
			float b1  = std_PrcpDiff_OBS/std_PrcpDiff_CON;
			float ro1 = -0.1;
			float ro2 = -0.1;

			float std_TavgDiff_BC = a1*std_TavgDiff_FUT;
			float std_PrcpDiff_BC = b1*std_PrcpDiff_FUT;

			if(opt_ecdf == 1)
			{
				//equations need to be revised
				fBCTDiffData[iDiffIndex] = std_TavgDiff_BC/std_TavgDiff_FUT * fTavgDiff 
											+ std_TavgDiff_FUT/std_PrcpDiff_FUT*ro1*fPrcpDiff;
				fBCPDiffData[iDiffIndex] = std_PrcpDiff_BC/std_PrcpDiff_FUT * fPrcpDiff
											+ std_PrcpDiff_FUT/std_TavgDiff_FUT*ro2*fTavgDiff;
			}
			else
			{
				fBCTDiffData[iDiffIndex] = fTavgDiff/std_TavgDiff_FUT * std_TavgDiff_BC;
				fBCPDiffData[iDiffIndex] = fPrcpDiff/std_PrcpDiff_FUT * std_PrcpDiff_BC;
			}
		}

	}//for(iDiffIndex)

	//////////////////////////////////////////////////////////////////
	//Step5: Reconsturcting Tavg and Prcp from
	//       the bias corrected differces sequences
	//////////////////////////////////////////////////////////////////
	//calculate climatology of obs/con/fut
	float fClimat[DATA_COL_NUM][12];
	for(int iCol=0; iCol<DATA_COL_NUM; iCol++)
	{
		for(int iMonth=0; iMonth<12; iMonth++)
		{
			float fTotalSum = 0.0f;
			for(int iData=0+iMonth; iData<iDataLengthArray[iCol]; iData+=12)
			{
				fTotalSum += fDataInput[iCol][iData];
			}
			fClimat[iCol][iMonth] = fTotalSum / (float) (iDataLengthArray[iCol]/12);
		}
	}

	//adding bias into pure incremental series
	for(int iData=0; iData<iDataLengthArray[4]; iData++)
	{
		fBCTData[iData] = fClimat[0][iData%12] + (fClimat[4][iData%12] - fClimat[2][iData%12]) + fBCTDiffData[iData];
		fBCPData[iData] = fClimat[1][iData%12] + (fClimat[5][iData%12] - fClimat[3][iData%12]) + fBCPDiffData[iData];
	}

	//make sure the precipitation is positive
	for(int i=0; i < iDataLengthArray[4]; i++)
	{
		if(fBCPData[i] < 0) 
		{
			//negtive extreme value processing
			fBCPData[i] = 1/(1-fBCPData[i]);
		}
	}

	////////////////////////////////////////////////////////////////////
	////Step6: Writting BC Sequences
	////////////////////////////////////////////////////////////////////
	for(int i=0; i < iDataLengthArray[4]; i++)
	{
		* (fBCTSTavg+i) = fBCTData[i];
		* (fBCTSPrcp+i) = fBCPData[i];
	}

	return 0;
}
