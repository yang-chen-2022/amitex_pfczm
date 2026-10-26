#!/bin/sh

comp_prefix=@CMAKE_INSTALL_PREFIX@
moddir=$comp_prefix/include
libdir=$comp_prefix/lib
version=@PROJECT_VERSION@

usage()
{
    cat <<EOF
Usage: FoX-config [OPTION]

Known values for OPTION are:

  --libs		print library linking information
  --fcflags		print pre-processor and compiler flags
  --help		display this help and exit
  --version		display version information
EOF

    exit $1
}

FLAGS=both
MODULES=
while test $# -gt 0; do
    case "$1" in
    -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) optarg= ;;
    esac

    case "$1" in
    --version)
	echo ${version}
	exit 0
	;;

    --help)
	usage 0
	;;

    --fcflags)
        FLAGS=FC
       	;;

    --libs)
        FLAGS=LD
       	;;

    --common)
        MODULES="$MODULES common utils fsys"
        ;;
    --utils)
        MODULES="$MODULES utils fsys"
        ;;
    --wxml)
        MODULES="$MODULES wxml common utils fsys"
        ;;
    --wcml)
        MODULES="$MODULES wcml wxml common utils fsys"
        ;;
    --sax)
        MODULES="$MODULES sax common utils fsys"
        ;;
    --dom)
        MODULES="$MODULES dom sax wxml common utils fsys"
        ;;
    --wkml)
        MODULES="$MODULES wkml wxml common utils fsys"
        ;;

    *)
	usage
	exit 1
	;;
    esac
    shift
done

if [ "x$MODULES" = x ]; then
#  MODULES="dom sax wcml wxml common utils fsys"
  MODULES="dom sax wkml wcml wxml common utils fsys"
fi

case $MODULES in *fsys*)
  LIBS="-lFoX_fsys"
  ;;
esac
case $MODULES in *utils*)
  LIBS="-lFoX_utils $LIBS"
  ;;
esac
case $MODULES in *common*)
  LIBS="-lFoX_common $LIBS"
  ;;
esac
case $MODULES in *wxml*)
  LIBS="-lFoX_wxml $LIBS"
  ;;
esac
case $MODULES in *wcml*)
  LIBS="-lFoX_wcml $LIBS"
  ;;
esac
case $MODULES in *sax*)
  LIBS="-lFoX_sax $LIBS"
  ;;
esac
case $MODULES in *dom*)
  LIBS="-lFoX_dom $LIBS"
  ;;
esac
case $MODULES in *wkml*)
  LIBS="-lFoX_wkml $LIBS"
  ;;
esac


case $FLAGS in
  FC)
    echo @FC_MOD_FLAG@"$moddir"
    ;;
  LD)
    echo $LIBS
    ;;
  both)
    echo @FC_MOD_FLAG@"$moddir" $LIBS
esac 

exit 0
