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

typedef void (*UMATPtr)(UMATReal *const,
			UMATReal *const,
			UMATReal *const,
			const UMATReal *const,
			const UMATReal *const,
			const UMATReal *const,
			const UMATReal *const,
			const UMATReal *const,
			const UMATReal *const,
			const UMATInt  *const,
			const UMATReal *const,
			const UMATInt  *const,
			UMATInt  *const);


void umatD_call(void * f,
	       UMATReal *const FLUXD,
	       UMATReal *const GRADQD,
	       UMATReal *const DGRADQD,
	       const UMATReal *const TIME,
	       const UMATReal *const DTIME,
	       const UMATReal *const TEMP,
	       const UMATReal *const DTEMP,
	       const UMATReal *const PREDEF,
	       const UMATReal *const DPRED,
	       const UMATInt  *const NVARD,
	       const UMATReal *const PROPS,
	       const UMATInt  *const NPROPS,
	       UMATInt *const KINC)
{
  union{
    void * f;
    UMATPtr ptr;
  } ptr;
  ptr.f = f;
  UMATPtr u=ptr.ptr;
  u(FLUXD,GRADQD,DGRADQD,TIME,DTIME,TEMP,DTEMP,  
       PREDEF,DPRED,NVARD,PROPS,NPROPS,KINC);
}
