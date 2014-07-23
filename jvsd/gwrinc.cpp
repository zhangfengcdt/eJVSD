// gwrinc.c
// functions of regird algorithm
// Feng Zhang <fzhang@gatech.edu>
// Mar., 2009

#include "stdafx.h"
#include "gwrinc.h"
#include <stdlib.h>
#include <stdio.h>


////////////////////////////
//Retirve Lon&Lat Dim Info
////////////////////////////
int get_grid_coord_dim(char* charFileName, char* charLatName, char* charLonName, int* pDimLat, int* pDimLon)
{
	if(!charFileName || !charLatName || !charLonName || !pDimLat || !pDimLon)
	{
		printf("Error: load_grid_coord.\n");
		return -1;
	}

	int ncid;
	int lat_varid, lon_varid;
	size_t latlength, lonlength;

	/* Error handling. */
	int retval;

	/* Open the file. */
	if ((retval = nc_open(charFileName, NC_NOWRITE, &ncid)))
		ERR(retval);

	/* Get the dimids of the latitude and longitude coordinate
	* variables. */
	if ((retval = nc_inq_dimid(ncid, charLatName, &lat_varid)))
		ERR(retval);
	if ((retval = nc_inq_dimid(ncid, charLonName, &lon_varid)))
		ERR(retval);

	/* Read the coordinate dimension length. */
	if ((retval = nc_inq_dimlen(ncid, lat_varid, &latlength)))
		ERR(retval);
	if ((retval = nc_inq_dimlen(ncid, lon_varid, &lonlength)))
		ERR(retval);

	*pDimLat = (int)latlength;
	*pDimLon = (int)lonlength;

	/* Close the file */
	if ((retval = nc_close(ncid)))
		ERR(retval);

	return 0;
}


////////////////////////////
//load grid coordinate
////////////////////////////
int load_grid_coord(char* charFileName, char* charLatName, char* charLonName, float* plat, float* plon)
{
	if(!charFileName || !charLatName || !charLonName || !plat || !plon)
	{
		printf("Error: load_grid_coord.\n");
		return -1;
	}

	int ncid;
	int lat_varid, lon_varid;

	/* Error handling. */
	int retval;

	/* Open the file. */
	if ((retval = nc_open(charFileName, NC_NOWRITE, &ncid)))
		ERR(retval);

	/* Get the varids of the latitude and longitude coordinate
	* variables. */
	if ((retval = nc_inq_varid(ncid, charLatName, &lat_varid)))
		ERR(retval);
	if ((retval = nc_inq_varid(ncid, charLonName, &lon_varid)))
		ERR(retval);

	/* Read the coordinate variable data. */
	if ((retval = nc_get_var_float(ncid, lat_varid, plat)))
		ERR(retval);
	if ((retval = nc_get_var_float(ncid, lon_varid, plon)))
		ERR(retval);

	/* Close the file */
	if ((retval = nc_close(ncid)))
		ERR(retval);

	return 0;
}

////////////////////////////
//load netCDF variable data
//for a given time slice
////////////////////////////
int load_netcdf_2dvar(char* charFileName, char* charVariableName, int time_id, int nLat, int nLon, float* pData)
{
	if(!charFileName || !charVariableName || !pData)
	{
		printf("Error: load_grid_coord.\n");
		return -1;
	}

	int ncid;
	int var_varid;

	/* Error handling. */
	int retval;

	/* Loop indexes. */
	int i = 0;

	/* The start and count arrays will tell the netCDF library where to
	read our data. */
	size_t start[3], count[3];

	/* Open the file. */
	if ((retval = nc_open(charFileName, NC_NOWRITE, &ncid)))
		ERR(retval);

	/* Get the varids of the netCDF variable. */
	if ((retval = nc_inq_varid(ncid, charVariableName, &var_varid)))
		ERR(retval);

	/* Read the data. Since we know the contents of the file we know
	* that the data arrays in this program are the correct size to
	* hold one timestep. */
	count[0] = 1;
	count[1] = nLat;
	count[2] = nLon;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;

	/* Read and check one record at a time. */
	start[0] = time_id;
	if ((retval = nc_get_vara_float(ncid, var_varid, start, 
		count, pData)))
		ERR(retval);

	/* Close the file */
	if ((retval = nc_close(ncid)))
		ERR(retval);

	return 0;
}

////////////////////////////
//save netCDF variable data
//for a given time slice
////////////////////////////
int save_netcdf_2dvar(char* charFileName, char* charVariableName, int nLat, int nLon, 
					  float* pData, char *charVarUnit, float* plat, float* plon, float fMissingValue)
{
	if(!charFileName || !charVariableName || !plat || !plon || !pData )
	{
		printf("Error: save_netcdf_2dvar.\n");
		return -1;
	}

   /* IDs for the netCDF file, dimensions, and variables. */
   int ncid, lon_dimid, lat_dimid;
   int lat_varid, lon_varid, var_varid;
   int dimids[2];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[2], count[2];

   /* Loop indexes. */
   int i = 0;
   
   /* Error handling. */
   int retval;

   /* Create the file. */
   if ((retval = nc_create(charFileName, NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions.*/
   if ((retval = nc_def_dim(ncid, LAT_NAME, nLat, &lat_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, LON_NAME, nLon, &lon_dimid)))
      ERR(retval);

   /* Define the coordinate variables. */
   if ((retval = nc_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid, 
			    &lat_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid, 
			    &lon_varid)))
      ERR(retval);

   /* Assign units attributes to coordinate variables. */
   if ((retval = nc_put_att_text(ncid, lat_varid, COOR_UNITS, 
				 strlen(DEGREES_NORTH), DEGREES_NORTH)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, lon_varid, COOR_UNITS, 
				 strlen(DEGREES_EAST), DEGREES_EAST)))
      ERR(retval);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = lat_dimid;
   dimids[1] = lon_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, charVariableName, NC_FLOAT, 2, 
			    dimids, &var_varid)))
      ERR(retval);

   /* Assign units attributes to the netCDF variables. */
   float fAttrValue = 0.0f;
   if ((retval = nc_put_att_text(ncid, var_varid, COOR_UNITS, 
				 strlen(charVarUnit), charVarUnit)))
      ERR(retval);

   if ((retval = nc_put_att_float(ncid, var_varid, "_FillValue", 
				 NC_FLOAT, 1, &fMissingValue)))
      ERR(retval);

   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_float(ncid, lat_varid, plat)))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid, lon_varid, plon)))
      ERR(retval);

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = nLat;
   count[1] = nLon;
   start[0] = 0;
   start[1] = 0;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
   if ((retval = nc_put_vara_float(ncid, var_varid, start, count, 
				      pData)))
   ERR(retval);

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

	return 0;
}

////////////////////////////
//load time coordinate
////////////////////////////
int load_time_coord(char* charFileName, char* charTimeName, int* ptime)
{
	if(!charFileName || !charTimeName || !ptime)
	{
		printf("Error: load_time_coord.\n");
		return -1;
	}

	int ncid;
	int time_varid;

	/* Error handling. */
	int retval;

	/* Open the file. */
	if ((retval = nc_open(charFileName, NC_NOWRITE, &ncid)))
		ERR(retval);

	/* Get the varids of the latitude and longitude coordinate
	* variables. */
	if ((retval = nc_inq_varid(ncid, charTimeName, &time_varid)))
		ERR(retval);

	/* Read the coordinate variable data. */
	if ((retval = nc_get_var_int(ncid, time_varid, ptime)))
		ERR(retval);

	/* Close the file */
	if ((retval = nc_close(ncid)))
		ERR(retval);

	return 0;
}


////////////////////////////
//create netCFD file (3D)
////////////////////////////
int create_netcdf_3dvar(char* charFileName, char* charVariableName, int nLat, int nLon, int nTime,
					  float* pData, char *charVarUnit, float* plat, float* plon, int* ptime, 
					  float fMissingValue, int *pncid, int *pvar_varid)
{
	if(!charFileName || !charVariableName || !charVarUnit || !plat || !plon 
		             || !pData || !pncid || !pvar_varid || !ptime)
	{
		printf("Error: create_netcdf_3dvar.\n");
		return -1;
	}

   /* IDs for the netCDF file, dimensions, and variables. */
   int ncid, lon_dimid, lat_dimid, time_dimid;
   int lat_varid, lon_varid, time_varid, var_varid;
   int dimids[3];

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
//   size_t start[3], count[3];

   /* Loop indexes. */
   int i = 0;
   
   /* Error handling. */
   int retval;

   /* Create the file. */
   if ((retval = nc_create(charFileName, NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions.*/
   if ((retval = nc_def_dim(ncid, REC_NAME, nTime,&time_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, LAT_NAME, nLat, &lat_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, LON_NAME, nLon, &lon_dimid)))
      ERR(retval);

   /* Define the coordinate variables. */
   if ((retval = nc_def_var(ncid, REC_NAME, NC_INT,  1, &time_dimid, 
			    &time_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid, 
			    &lat_varid)))
      ERR(retval);
   if ((retval = nc_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid, 
			    &lon_varid)))
      ERR(retval);

   /* Assign units attributes to coordinate variables. */
   if ((retval = nc_put_att_text(ncid, time_varid, TIME_UNITS, 
				 strlen(TIME_CRU), TIME_CRU)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, lat_varid, COOR_UNITS, 
				 strlen(DEGREES_NORTH), DEGREES_NORTH)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid, lon_varid, COOR_UNITS, 
				 strlen(DEGREES_EAST), DEGREES_EAST)))
      ERR(retval);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = time_dimid;
   dimids[1] = lat_dimid;
   dimids[2] = lon_dimid;

   /* Define the netCDF variables for the pressure and temperature
    * data. */
   if ((retval = nc_def_var(ncid, charVariableName, NC_FLOAT, 3, 
			    dimids, &var_varid)))
      ERR(retval);

   /* Assign units attributes to the netCDF variables. */
   if ((retval = nc_put_att_text(ncid, var_varid, COOR_UNITS, 
				 strlen(charVarUnit), charVarUnit)))
      ERR(retval);

   if ((retval = nc_put_att_float(ncid, var_varid, "_FillValue", 
				 NC_FLOAT, 1, &fMissingValue)))
      ERR(retval);

   /* End define mode. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_int(ncid, time_varid, ptime)))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid, lat_varid, plat)))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid, lon_varid, plon)))
      ERR(retval);

   *pncid       = ncid;
   *pvar_varid  = var_varid;

   return 0;
}

////////////////////////////
//save netCDF variable data
////////////////////////////
int save_netcdf_3dvar(float* pData, int nTime, int nLat, int nLon, int iTime, int ncid, int var_varid)
{
	if(!pData)
	{
		printf("Error: save_netcdf_3dvar.\n");
		return -1;
	}

   /* Error handling. */
   int retval;

   /* The start and count arrays will tell the netCDF library where to
      write our data. */
   size_t start[3], count[3];

   /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
   count[0] = 1;
   count[1] = nLat;
   count[2] = nLon;
   start[0] = 0;
   start[1] = 0;
   start[2] = 0;

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays only hold one timestep worth
      of data. We will just rewrite the same data for each timestep. In
      a real application, the data would change between timesteps. */
   start[0] = iTime;
   if ((retval = nc_put_vara_float(ncid, var_varid, start, count, 
				      pData)))
   ERR(retval);
   return 0;
}


////////////////////////////
//close netCDF variable data
////////////////////////////
int close_netcdf_3dvar(int ncid)
{
   /* Error handling. */
   int retval;

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);
	
   return 0;
}
