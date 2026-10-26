#include<stdio.h>
#include<string.h>
#include<stdlib.h>

// Définition des types
//#ifdef UNIX32
//typedef int   UMATInt;
//#endif

//#ifdef UNIX64
//#ifdef WIN64
//typedef long long UMATInt;
//#else
//typedef long UMATInt;
//#endif
//#endif

// modif LG 7/8/18: 
typedef int64_t UMATInt;

typedef double UMATReal;

typedef void (*UMATPtr)(UMATReal *const,UMATReal *const,
			UMATReal *const,UMATReal *const,
			UMATReal *const,UMATReal *const,
			UMATReal *const,UMATReal *const,
			UMATReal *const,UMATReal *const,
			const UMATReal *const,const UMATReal *const,
			const UMATReal *const,const UMATReal *const,
			const UMATReal *const,const UMATReal *const,
			const UMATReal *const,const UMATReal *const,
			const char     *const,const UMATInt  *const,
			const UMATInt  *const,const UMATInt  *const,
			const UMATInt  *const,const UMATReal *const,
			const UMATInt  *const,const UMATReal *const,
			const UMATReal *const,const UMATReal *const,
			const UMATReal *const,const UMATReal *const,
			const UMATReal *const,const UMATInt  *const,
			const UMATInt  *const,const UMATInt  *const,
			const UMATInt  *const,const UMATInt  *const,
			UMATInt  *const,const int);


void umat_call(void * f,
	       UMATReal *const STRESS,
	       UMATReal *const STATEV,
	       UMATReal *const DDSDDE,
	       UMATReal *const SSE,
	       UMATReal *const SPD,
	       UMATReal *const SCD,
	       UMATReal *const RPL,
	       UMATReal *const DDSDDT,
	       UMATReal *const DRPLDE,
	       UMATReal *const DRPLDT,
	       const UMATReal *const STRAN,
	       const UMATReal *const DSTRAN,
	       const UMATReal *const TIME,
	       const UMATReal *const DTIME,
	       const UMATReal *const TEMP,
	       const UMATReal *const DTEMP,
	       const UMATReal *const PREDEF,
	       const UMATReal *const DPRED,
	       const char     *const CMNAME,
	       const UMATInt  *const NDI,
	       const UMATInt  *const NSHR,
	       const UMATInt  *const NTENS,
	       const UMATInt  *const NSTATV,
	       const UMATReal *const PROPS,
	       const UMATInt  *const NPROPS,
	       const UMATReal *const COORDS,
	       const UMATReal *const DROT,
	       const UMATReal *const PNEWDT,
	       const UMATReal *const CELENT,
	       const UMATReal *const DFGRD0,
	       const UMATReal *const DFGRD1,
	       const UMATInt  *const NOEL,
	       const UMATInt  *const NPT,
	       const UMATInt  *const LAYER,
	       const UMATInt  *const KSPT,
	       const UMATInt  *const KSTEP,
	       UMATInt *const KINC)
{
  union{
    void * f;
    UMATPtr ptr;
  } ptr;
  ptr.f = f;
  UMATPtr u=ptr.ptr;
  u(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,
    DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
    CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
    PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,
    KINC,strlen(CMNAME));
}
