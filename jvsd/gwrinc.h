// gwrinc.h
// functions of NetCDF operations
// Feng Zhang <fzhang@gatech.edu>
// Mar., 2009

#ifndef _GWRI_NC_
#define _GWRI_NC_

#include <stddef.h> /* size_t, ptrdiff_t */
#include <errno.h>  /* netcdf functions sometimes return system errors */
#include <netcdf.h>

#if defined(__cplusplus)
extern "C" {
#endif

/* Handle errors by printing an error message and exiting with a
* non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

/*
 *  The netcdf external data types
 */
//typedef enum {
//	REGRID_ERROR_NONUNI     =	0,
//	REGRID_ERROR_NULLBUF    =	1,
//	REGRID_ERROR_BUFLENG    =   2 
//} sample_enum;

/* For dimension names. */
#define REC_NAME "time"
#define LAT_NAME "lat"
#define LON_NAME "lon"

/* For the units attributes. */
#define COOR_UNITS    "units"
#define TIME_UNITS    "units"
#define DEGREES_EAST  "degrees_east"
#define DEGREES_NORTH "degrees_north"
#define TIME_CRU      "months since 1900-01-15 00:00:00"

/*
 * The Interface
 */
# define EXTERNL_ extern

//EXTERNL_ int
//gwri_stat_climatology(float* high_lat, float* high_lon, int high_lat_dim, int high_lon_dim, 
//						    float* low_lat,  float* low_lon,  int low_lat_dim, int low_lon_dim, 
//						    float* pin, float* pout, float fMissingValue);


////////////////////////////
//Retirve Lon&Lat Dim Info
////////////////////////////
EXTERNL_ int
get_grid_coord_dim(char* charFileName, char* charLatName, char* charLonName, int* pDimLat, int* pDimLon);

////////////////////////////
//load grid coordinate
////////////////////////////
EXTERNL_ int
load_grid_coord(char* charFileName, char* charLatName, char* charLonName, float* plat, float* plon);

////////////////////////////
//load netCDF variable data
//for a given time slice
////////////////////////////
EXTERNL_ int
load_netcdf_2dvar(char* charFileName, char* charVariableName, int time_id, int nLat, int nLon, float* pData);

////////////////////////////
//save netCDF variable data
//for a given time slice
////////////////////////////
EXTERNL_ int
save_netcdf_2dvar(char* charFileName, char* charVariableName, int nLat, int nLon, 
					  float* pData, char *charVarUnit, float* plat, float* plon, float fMissingValue);

////////////////////////////
//load time coordinate
////////////////////////////
EXTERNL_ int
load_time_coord(char* charFileName, char* charTimeName, int* ptime);


////////////////////////////
//create netCFD file (3D)
////////////////////////////
EXTERNL_ int
create_netcdf_3dvar(char* charFileName, char* charVariableName, int nLat, int nLon, int nTime,
					  float* pData, char *charVarUnit, float* plat, float* plon, int* ptime, 
					  float fMissingValue, int *pncid, int *pvar_varid);

////////////////////////////
//save netCDF variable data
////////////////////////////
EXTERNL_ int
save_netcdf_3dvar(float* pData, int nTime, int nLat, int nLon, int iTime, int ncid, int var_varid);


////////////////////////////
//close netCDF variable data
////////////////////////////
EXTERNL_ int
close_netcdf_3dvar(int ncid);


#if defined(__cplusplus)
}
#endif

#endif /* _GWRI_NC_ */
