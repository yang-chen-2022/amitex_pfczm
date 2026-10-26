#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <dlfcn.h>

void * umat_load(const char * const l,
		  const char * const n)
{
  void *lib;
  void *fct;
  lib = dlopen(l,RTLD_NOW);
  if(lib==0){
    fprintf(stderr,"umat_load: loading library '%s'"
	    " failed (%s)\n",l,dlerror());
    exit(-1);
  }
  fct = dlsym(lib,n);
  if(fct==0){
    fprintf(stderr,"umat_load: loading function '%s'"
	    " failed (%s)\n",n,dlerror());
    exit(-1);
  }
//  fprintf(stdout,"fct: %p\n",fct);
  return fct;
}
