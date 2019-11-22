#!/bin/bash
CPO=$JULIA_DEPOT_PATH/cpo
GUROBI=$CPO/gurobi
IPOPT=$CPO/ipopt
INTEL=$CPO/intel
SCIPOPTDIR=$CPO/scipoptsuite
LICENSE=$CPO/license
mkdir -p $SCIPOPTDIR $LICENSE && cd /tmp

# gurobi
if [ ! -f $GUROBI/linux64/bin/grbgetkey ]; then
    wget https://packages.gurobi.com/8.0/gurobi8.0.1_linux64.tar.gz
    tar xvzf gurobi*.tar.gz && mv gurobi801 $GUROBI && \rm gurobi*.tar.gz
    $GUROBI/*/*/grbgetkey -q --path $LICENSE 5b0babc0-0a23-11ea-9bcf-0a7c4f30bdbe
fi

# mkl
if [ ! -f $INTEL/bin/compilervars.sh ]; then
    wget http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/13575/l_mkl_2019.0.117.tgz
    tar xvzf l_mkl_*.tgz && cd l_mkl_2019.0.117
    sed -i "s:PSET_INSTALL_DIR=/opt/intel:PSET_INSTALL_DIR=$INTEL:g" silent.cfg
    sed -i "s:ARCH_SELECTED=ALL:ARCH_SELECTED=INTEL64:g" silent.cfg
    echo "ACCEPT_EULA=accept" >> silent.cfg
    ./install.sh --silent silent.cfg
    cd /tmp && \rm -rf l_mkl_*
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
    make && make install && cd /tmp && \rm -rf CoinIpopt
fi

# scip
if [ ! -f $SCIPOPTDIR/bin/scip ]; then
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
    make && make install && cd /tmp && \rm -rf scipoptsuite-*
fi

# license
cat > $LICENSE/mosek.lic <<- EOF
VENDOR MOSEKLM
FEATURE PTS MOSEKLM 9 17-nov-2020 uncounted VENDOR_STRING=-1*-1*-1*0 \
	HOSTID=DEMO TS_OK SIGN="1B99 AA0C B1EB 2317 A384 7390 187A \
	2241 8A56 1863 B70A A62A 51C2 14D1 E90E 0BDC BB79 115F 1E34 \
	8BDA 2BBA 32E6 9A5D 777B FE01 CDE8 356C 2541 5124 5F16" \
	SIGN2="1201 FA9A DBA0 0A88 53FB 34B7 4C80 ABEA CD19 39EB 5644 \
	1B7C 0200 7156 689D 186F 738F 9CDE CB6B C222 7EFA A60D 792D \
	A4F4 3B97 E6A0 C5F9 5394 BACD E269"
FEATURE PTON MOSEKLM 9 17-nov-2020 uncounted VENDOR_STRING=-1*-1*-1*0 \
	HOSTID=DEMO TS_OK SIGN="18A8 09C0 3CD2 09FA B68B 8877 A8DD \
	6F1B ECD4 AD4F 3E98 287D 95D2 DA3E 7B0C 058D D9BE F1D9 902E \
	04C8 DA6C 865E 7633 352B 6E14 EDE0 3BEB DA62 A7C1 3F60" \
	SIGN2="1684 BCF7 D5AE 904F 5E71 8C43 43DB 5AE9 4B02 56DF 44A4 \
	7624 BB2F D2CC D174 109D 1C7F 78CA DB92 4CFC 6E63 520E A712 \
	C900 7054 4075 A568 9478 EF1D E9E5"
EOF

cat > $LICENSE/pardiso.lic <<< '7B0D336BF36A2751AF3B391A642D903EFEDD038F1DE83F43FF6F909B'

cat > $LICENSE/wsmp.lic <<- EOF
-1045953726
-325865645
1736384656
-1939018359
2075595150
371191215
1784526268
581818635
1084814312
-1884243425
-638683284
1014231349
349833674
738256699
-713515432
1957325745
-716834666
2056043543
1610958852
893981165
-35583966
-1785265485
1032904816
1639046912
1152262969
66489726
831181535
EOF