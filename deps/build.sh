#!/bin/bash
CWD=$(pwd)
GUROBI=$CWD/usr/gurobi
IPOPT=$CWD/usr/ipopt
INTEL=$CWD/usr/intel
SCIPOPTDIR=$CWD/usr/scipoptsuite
FAKETIME=$CWD/usr/libfaketime
mkdir -p $SCIPOPTDIR

# gurobi
if [ ! -f $GUROBI/linux64/bin/grbgetkey ]; then
    wget https://packages.gurobi.com/8.0/gurobi8.0.1_linux64.tar.gz
    tar xvzf gurobi*.tar.gz && mv gurobi801 $GUROBI && \rm gurobi*.tar.gz
fi

# mkl
if [ ! -f $INTEL/bin/compilervars.sh ]; then
    wget http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/13575/l_mkl_2019.0.117.tgz
    tar xvzf l_mkl_*.tgz && cd l_mkl_2019.0.117
    sed -i "s:PSET_INSTALL_DIR=/opt/intel:PSET_INSTALL_DIR=$INTEL:g" silent.cfg
    sed -i "s:ARCH_SELECTED=ALL:ARCH_SELECTED=INTEL64:g" silent.cfg
    echo "ACCEPT_EULA=accept" >> silent.cfg
    ./install.sh --silent silent.cfg
    cd $CWD && \rm -rf l_mkl_*
fi
source $INTEL/bin/compilervars.sh intel64

# ipopt
if [ ! -f $IPOPT/bin/ipopt ]; then
    svn co https://projects.coin-or.org/svn/Ipopt/stable/3.12 CoinIpopt
    cd CoinIpopt/ThirdParty
    cd HSL && wget https://sourceforge.net/projects/bearapps/files/coinhsl-2015.06.23.tar.gz \
    && tar xvzf coinhsl*.tar.gz && rm coinhsl*.tar.gz && mv coinhsl-2015.06.23 coinhsl && cd ..
    cd ASL && ./get.ASL && cd ..
    cd Blas && ./get.Blas && cd ..
    cd Lapack && ./get.Lapack && cd ..
    cd Metis && ./get.Metis && cd ..
    cd Mumps && ./get.Mumps && cd ..
    cd .. && mkdir -p build && cd build
    ../configure --prefix=$IPOPT coin_skip_warn_cxxflags=yes \
        ADD_FFLAGS=-fPIC ADD_CFLAGS=-fPIC ADD_CXXFLAGS=-fPIC \
        ADD_CFLAGS=-fopenmp ADD_FFLAGS=-fopenmp ADD_CXXFLAGS=-fopenmp
        --with-blas="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl"
    make && make install && cd $CWD && \rm -rf CoinIpopt
fi

# scip
if [ -f $SCIPOPTDIR/bin/scip ]; then
    wget https://download.mosek.com/stable/8.1.0.64/mosektoolslinux64x86.tar.bz2
    tar xvjf mosek*.tar.bz2 && \rm mosek*.tar.bz2 && mv mosek $SCIPOPTDIR
    MSKDIR=$(echo $SCIPOPTDIR/mosek/*/tools/platform/*)
    wget https://sourceforge.net/projects/bearapps/files/scipoptsuite-6.0.0.tgz
    tar xvzf scipoptsuite-*.tgz && \rm scipoptsuite-*.tgz && cd scipoptsuite-*
    wget https://sourceforge.net/projects/bearapps/files/worhp_1.12-3_linux.zip
    unzip -o -d $SCIPOPTDIR/worhp worhp_1.12-3_linux.zip
    mkdir -p build && cd build
    export CPATH=$CPATH:$IPOPT/include/coin
    cmake -DCMAKE_INSTALL_PREFIX=$SCIPOPTDIR -DGMP=off -DREADLINE=off -DPARASCIP=on \
        -DLPS=msk -DMOSEK_INCLUDE_DIRS=$MSKDIR/h -DMOSEK_LIBRARY=$MSKDIR/bin/libmosek64.so \
        -DWORHP=off -DWORHP_INCLUDE_DIRS=$SCIPOPTDIR/worhp/include -DWORHP_LIBRARY=$SCIPOPTDIR/worhp/lib/libworhp.so \
        -DIPOPT=on -DIPOPT_LIBRARIES=$IPOPT/lib/libipopt.so ..
    make && make install && cd $CWD && \rm -rf scipoptsuite-*
fi

# faketime
if [ -f $FAKETIME/src/libfaketime.so.1 ]; then
    git clone https://github.com/wolfcw/libfaketime.git $FAKETIME
    cd $FAKETIME/src && make
fi