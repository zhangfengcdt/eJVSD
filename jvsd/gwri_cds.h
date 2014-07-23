#define DS_RATIO_TEMP        1
#define DS_RATIO_PRCP        800
//#define DS_RATIO_PRCP        400
//#define DS_RATIO_PRCP        200
//#define DS_RATIO_PRCP        100
//#define DS_RATIO_PRCP        50
//#define DS_RATIO_PRCP        10

///////////////////////////////////
//Retrieve Overlapping cells's ID
///////////////////////////////////
int get_overlap_cells(float *fpHighResuLat,  int nHighResuLat,
					  float *fpHighResuLon,  int nHighResuLon,
					  float *fpLowResuLat,   int nLowResuLat,
					  float *fpLowResuLon,   int nLowResuLon,
					  float fLat,            float fLon,
					  int *npMaxOverlapLat, int *npMinOverlapLat,
					  int *npMaxOverlapLon, int *npMinOverlapLon,
					  float fExpandCoeff=0.5)
{
	if(!fpHighResuLat    || !fpHighResuLon   ||
	   !fpLowResuLat     || !fpLowResuLon    ||
	   !npMaxOverlapLat  || !npMinOverlapLat ||
	   !npMaxOverlapLat  || !npMinOverlapLon)
	{
		printf("Fatal error: (get_overlap_cells) pointer to TS is null.\n");
		return -1;
	}
	if(nHighResuLat < 0 || nHighResuLon < 0 ||
	   nLowResuLat < 0 || nLowResuLon < 0 )
	{
		printf("Fatal error: (ds_search) length of TS is negtive.\n");
		return -1;
	}

	//calculate the interval 
	//of high and low grid
	float high_lat_interv = 0.0f;
	float high_lon_interv = 0.0f;
	float low_lat_interv  = 0.0f;
	float low_lon_interv  = 0.0f;
	high_lat_interv = *(fpHighResuLat+10) - *(fpHighResuLat+9);
	high_lon_interv = *(fpHighResuLon+10) - *(fpHighResuLon+9);
	low_lat_interv  = *(fpLowResuLat+10)  - *(fpLowResuLat+9);
	low_lon_interv  = *(fpLowResuLon+10)  - *(fpLowResuLon+9);

	//get min,max of both grid
	float fHighLowerLat   = *fpHighResuLat;
	float fHighLowerLon   = *fpHighResuLon;
	if(fHighLowerLon<0) fHighLowerLon = 360.0f + fHighLowerLon;

	//calculate min,max of overlap
	//cells on high reslu. grid
	float fLeftCornerLat  = fLat - low_lat_interv*fExpandCoeff;
	float fLeftCornerLon  = fLon - low_lon_interv*fExpandCoeff;
	float fRightCornerLat = fLat + low_lat_interv*fExpandCoeff;
	float fRightCornerLon = fLon + low_lon_interv*fExpandCoeff;

	int iIntLook_LatMin = int((fLeftCornerLat  - fHighLowerLat) / high_lat_interv);
	int iIntLook_LatMax = int((fRightCornerLat - fHighLowerLat) / high_lat_interv);
	int iIntLook_LonMin = int((fLeftCornerLon  - fHighLowerLon) / high_lon_interv);
	int iIntLook_LonMax = int((fRightCornerLon - fHighLowerLon) / high_lon_interv);
	if(iIntLook_LatMin < 0) iIntLook_LatMin = 0;
	if(iIntLook_LatMax < 0) iIntLook_LatMax = 0;
	if(iIntLook_LonMin < 0) iIntLook_LonMin = 0;
	if(iIntLook_LonMax < 0) iIntLook_LonMax = 0;
	if(iIntLook_LatMin > nHighResuLat-1) iIntLook_LatMin = nHighResuLat-1;
	if(iIntLook_LatMax > nHighResuLat-1) iIntLook_LatMax = nHighResuLat-1;
	if(iIntLook_LonMin > nHighResuLon-1) iIntLook_LonMin = nHighResuLon-1;
	if(iIntLook_LonMax > nHighResuLon-1) iIntLook_LonMax = nHighResuLon-1;

	//return the values to caller
	*npMaxOverlapLat = iIntLook_LatMax;
	*npMinOverlapLat = iIntLook_LatMin;
	*npMaxOverlapLon = iIntLook_LonMax;
	*npMinOverlapLon = iIntLook_LonMin;

	return 1;
}

///////////////////////////////////
//Downscaling searching and return 
//the index of monthly sequences
//opt_flag = 0: P+T can be diff month
//opt_flag = 1: P+T must be same month
//Monthly version
///////////////////////////////////
int ds_search_monthly(float* fOBSTSTavg,      int nOBSTSTavg,
			  float* fOBSTSPrcp,      int nOBSTSPrcp,
			  float* fBCTSTavg,       int nBCTSTavg,
			  float* fBCTSPrcp,       int nBCTSPrcp,
			  float* fIndexDSTSTavg,  int nDSTSTavg,
			  float* fIndexDSTSPrcp,  int nDSTSPrcp,
			  int opt_flag, float fMissingValue)
{
	if(!fOBSTSTavg || !fOBSTSPrcp ||
	   !fIndexDSTSTavg || !fIndexDSTSPrcp ||
	   !fBCTSTavg || !fBCTSPrcp)
	{
		printf("Fatal error: (ds_search_monthly) pointer to TS is null.\n");
		return -1;
	}

	if(nOBSTSTavg < 0 || nOBSTSTavg < 0 ||
	   nDSTSTavg < 0 || nDSTSPrcp < 0 ||
	   nBCTSTavg  < 0 || nBCTSPrcp < 0)
	{
		printf("Fatal error: (ds_search_monthly) length of TS is negtive.\n");
		return -1;
	}
	if(nOBSTSTavg != nOBSTSPrcp || 
	   nDSTSTavg  != nDSTSPrcp||
	   nBCTSTavg  != nBCTSPrcp)
	{
		printf("Fatal error: (ds_search_monthly) lengthes of TSs are not equal.\n");
		return -1;
	}

	//Step1: find the lower and upper
	//bound of observational sequences
	float fMinOBSTavgMonthly[12], fMaxOBSTavgMonthly[12];
	float fMinOBSPrcpMonthly[12], fMaxOBSPrcpMonthly[12];
	float fMinOBSTavg=1e10, fMaxOBSTavg=-1e10;
	float fMinOBSPrcp=1e10, fMaxOBSPrcp=-1e10;
	for(int iMonth=0; iMonth<12; iMonth++)
	{
		for(int i=0; i<nOBSTSTavg; i++)
		{
			if(i%12 != iMonth) continue;
			if(*(fOBSTSTavg+i) < fMinOBSTavg) fMinOBSTavg = *(fOBSTSTavg+i);
			if(*(fOBSTSTavg+i) > fMaxOBSTavg) fMaxOBSTavg = *(fOBSTSTavg+i);
			if(*(fOBSTSPrcp+i) < fMinOBSPrcp) fMinOBSPrcp = *(fOBSTSPrcp+i);
			if(*(fOBSTSPrcp+i) > fMaxOBSPrcp) fMaxOBSPrcp = *(fOBSTSPrcp+i);
		}
		fMinOBSTavgMonthly[iMonth] = fMinOBSTavg;
		fMaxOBSTavgMonthly[iMonth] = fMaxOBSTavg;
		fMinOBSPrcpMonthly[iMonth] = fMinOBSPrcp;
		fMaxOBSPrcpMonthly[iMonth] = fMaxOBSPrcp;
	}

	//Step2: Look over each BC time sequences
	for(int i=0; i<nBCTSTavg; i++)
	{
		//Find the Month of the data
		int iBCMonth = i%12;

		if( *(fBCTSTavg+i) < fMinOBSTavgMonthly[iBCMonth] || *(fBCTSTavg+i) > fMaxOBSTavgMonthly[iBCMonth] ||
		    *(fBCTSPrcp+i) < fMinOBSPrcpMonthly[iBCMonth] || *(fBCTSPrcp+i) > fMaxOBSPrcpMonthly[iBCMonth] )
		{
			*(fIndexDSTSTavg + i) = fMissingValue;
			*(fIndexDSTSPrcp + i) = fMissingValue;
			continue;
		}

		float fbcT = *(fBCTSTavg+i), fbcP = *(fBCTSPrcp+i);
		float Rmin_square = 1e10;
		int   indexRmin   = 0;
		for(int j=0; j<nOBSTSTavg; j++)
		{
			//only look into the same month
			int iOBSMonth = j%12;
			if(iBCMonth != iOBSMonth) continue;

			float fobsT = *(fOBSTSTavg+j);
			float fobsP = *(fOBSTSPrcp+j);
			
			float R2 = DS_RATIO_TEMP*(fbcT-fobsT)*(fbcT-fobsT)+DS_RATIO_PRCP*(fbcP-fobsP)*(fbcP-fobsP);
			if(R2 < Rmin_square)
			{
				Rmin_square = R2;
				indexRmin   = j;
			}
		}

		*(fIndexDSTSTavg + i) = (float)indexRmin;
		*(fIndexDSTSPrcp + i) = (float)indexRmin;
	}

	return 1;
}

///////////////////////////////////
//Downscaling searching and return 
//the index of monthly sequences
//opt_flag = 0: P+T can be diff month
//opt_flag = 1: P+T must be same month
//check_bound = 0: do not check for outof bound value
//check_bound = 1: check for outof bound value
//Daily version
///////////////////////////////////
int ds_search_daily(float* fOBSTSTavg,      int nOBSTSTavg,
			  float* fOBSTSPrcp,      int nOBSTSPrcp,
			  float* fBCTSTavg,       int nBCTSTavg,
			  float* fBCTSPrcp,       int nBCTSPrcp,
			  float* fIndexDSTSTavg,  int nDSTSTavg,
			  float* fIndexDSTSPrcp,  int nDSTSPrcp,
			  int opt_flag, float fMissingValue, int check_bound)
{
	if(!fOBSTSTavg || !fOBSTSPrcp ||
	   !fIndexDSTSTavg || !fIndexDSTSPrcp ||
	   !fBCTSTavg || !fBCTSPrcp)
	{
		printf("Fatal error: (ds_search_daily) pointer to TS is null.\n");
		return -1;
	}

	if(nOBSTSTavg < 0 || nOBSTSTavg < 0 ||
	   nDSTSTavg < 0 || nDSTSPrcp < 0 ||
	   nBCTSTavg  < 0 || nBCTSPrcp < 0)
	{
		printf("Fatal error: (ds_search_daily) length of TS is negtive.\n");
		return -1;
	}
	if(nOBSTSTavg != nOBSTSPrcp || 
	   nDSTSTavg  != nDSTSPrcp||
	   nBCTSTavg  != nBCTSPrcp)
	{
		printf("Fatal error: (ds_search_daily) lengthes of TSs are not equal.\n");
		return -1;
	}

	//Step1: find the lower and upper
	//bound of observational sequences
	float fMinOBSTavg=1e10, fMaxOBSTavg=-1e10;
	float fMinOBSPrcp=1e10, fMaxOBSPrcp=-1e10;
	for(int i=0; i<nOBSTSTavg; i++)
	{
		if(*(fOBSTSTavg+i) < fMinOBSTavg) fMinOBSTavg = *(fOBSTSTavg+i);
		if(*(fOBSTSTavg+i) > fMaxOBSTavg) fMaxOBSTavg = *(fOBSTSTavg+i);
		if(*(fOBSTSPrcp+i) < fMinOBSPrcp) fMinOBSPrcp = *(fOBSTSPrcp+i);
		if(*(fOBSTSPrcp+i) > fMaxOBSPrcp) fMaxOBSPrcp = *(fOBSTSPrcp+i);
	}

	//Step2: Look over each BC time sequences
	for(int i=0; i<nBCTSTavg; i++)
	{
		if(*(fIndexDSTSTavg + i) !=fMissingValue 
		&& *(fIndexDSTSPrcp + i) !=fMissingValue)
		{
			continue;
		}

		if(check_bound == 1)
		{

			if( *(fBCTSTavg+i) < fMinOBSTavg || *(fBCTSTavg+i) > fMaxOBSTavg ||
				*(fBCTSPrcp+i) < fMinOBSPrcp || *(fBCTSPrcp+i) > fMaxOBSPrcp )
			{
				*(fIndexDSTSTavg + i) = fMissingValue;
				*(fIndexDSTSPrcp + i) = fMissingValue;
				continue;
			}
		}

		float fbcT = *(fBCTSTavg+i), fbcP = *(fBCTSPrcp+i);
		float Rmin_square = 1e10;
		int   indexRmin   = 0;
		for(int j=0; j<nOBSTSTavg; j++)
		{
			float fobsT = *(fOBSTSTavg+j);
			float fobsP = *(fOBSTSPrcp+j);
			float R2 = DS_RATIO_TEMP*(fbcT-fobsT)*(fbcT-fobsT)+DS_RATIO_PRCP*(fbcP-fobsP)*(fbcP-fobsP);
			if(R2 < Rmin_square)
			{
				Rmin_square = R2;
				indexRmin   = j;
			}
		}

		//return negtive value of index 
		//for daily sequences
		*(fIndexDSTSTavg + i) = (float)(-indexRmin);
		*(fIndexDSTSPrcp + i) = (float)(-indexRmin);
	}

	return 1;
}

