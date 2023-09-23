#!/bin/bash
##!/bin/bash -x
##!/bin/bash -xvD
#===========================================================
#
#    file:       prepare_case.sh
#
#    purpose:    bash script used to prepare (and potentially run) the BEPS model
#                for a new test/experiment.
#
#    created:    2020/03
#    last:       2022/07
#
#    author:     Michael Vossbeck (The Inversion Lab)
#
#=============================




dir_sitecmp_netcdf() {
  dir1=$1
  dir2=$2
  for ncf1 in `ls ${dir1}/beps_site_*.nc`
  do
    ncf2=${dir2}/`basename ${ncf1}`
    diff ${ncf1} ${ncf2} > /dev/null 2>&1
    if [ $? -ne 0 ]
    then
      echo "ERROR::DIFFERING FILES   ***${ncf1}***   ***${ncf2}***"
      exit 1
    fi
  done
}


##################################################
#
#--        d e f a u l t   s e t t i n g s
#
script_name=$0
basedir=$(readlink -f ${script_name} | xargs dirname | xargs dirname)
#-- relevant directories
script_dir=${basedir}/scripts
bld_dir=${basedir}/bld
model_dir=${basedir}/models_assim
ipt_dir=${basedir}/inputdata
#-- run time configuration (most depend on settings made on command line)
initial_yyyymmdd="20130101"
ndays=365
#>> override file name prefixes for lai/boundary/meteo forcing files
filename_prefixes=""
pert=0.
runc=-1
disable_prior=0
#--
external_nrlib=
#--
comparison_exe=${script_dir}/beps_inspect.py
exp_dir=
out_dir=
is_vanilla=0
tst_org=0
tst_runsimobs=0
tst_twin=0
tst_jacfwv=0
#-- whether to compile with debug flags enabled
dbg=0
#-- observational data file
assim_obs=



#---------------------------------------
#         usage
#
usage() {
    echo "------------------------------"
    echo
    echo "  usage: $0 [opts]"
    printf "  %-20s%-s\n" "--help" \
	   "print this message"
    printf "  %-20s%-s\n" "--vanilla"\
	   "prepare directory for compiling and running a new BEPS experiment."
    printf "  %-20s%-s\n" "--org"\
	   "create experiment directory, compile original BEPS driver, and run exectuable"
    printf "  %-20s%-s\n" "--runsimobs"\
	   "prepare experiment directory, compile BEPS driver that generates pseudo-observations, and run exectuable."
    printf "  %-20s%-s\n" "--twin"\
	   "prepare experiment directory, compile codes, and run an identical twin experiment."
    printf "  %-20s%-s\n" "--jacfwv"\
	   "compile derivative code driver, run BEPS vector-tangent-linear model and write Jacobian to ASCII and binary file."
    printf "  %-20s%-s\n" "--assim_obs"\
	   "provide an observational data file (NetCDF in consistent format) with fully qualified filepath."\
	   "If given, prepare experiment directory, compile minimisation driver, and run the data assimilation experiment. "\
	   "NOTE, the user is responsible to ensure consistency of observational spatio-temporal settings with the forcing files, "\
	   "and in particular ensure that settings made with options "\
	   "--yyyymmdd, --ndays, --nlp (or their default values) are consisten "\
	   "with the real observations provided!"\
	   "To make sure to start minimisation from prior control vector check that the perturbation (provided by --pert) has a value of zero."
    echo
    printf "  %-20s%-s\n" "--nrlib"\
	   "provide path to pre-built Numerical Recipes library (must be compliant with compiler envirionment with which the experiment will be compiled)"
    echo
    printf "  %-20s%-s\n" "--initial_yyyymmdd"\
	   "first day of simulation period (default: ${initial_yyyymmdd})"
    printf "  %-20s%-s\n" "--ndays"\
	   "length of simulation in days (default: ${ndays})"
    printf "  %-20s%-s\n" "--nlp" \
	   "number of sites or grid points (default:${nlp})."
    printf "  %-20s%-s\n" "--filename_prefixes" \
	   "provide file name prefixes for lai/boundary/meteo input files passed to the namelist builder"
    printf "  %-20s%-s\n" "--pert" \
	   "perturbation of initial control vector (for twin experiment) (default:${pert})."
    printf "  %-20s%-s\n" "--runc" \
	   "generate pseudo observations with relative uncertainty (negative values indicate uniform absolute uncertainty of 1, default:${runc})."
    printf "  %-20s%-s\n" "--disable_prior" \
	   "whether to disable prior-term in cost function (potentially for twin experiment) (default:${disable_prior})."
    echo
    printf "  %-20s%-s\n" "--expdir"\
	   "basic directory for the new BEPS experiment (suitable default will be created if none is given)"
    printf "  %-20s%-s\n" "--outdir"\
	   "output directory where BEPS generated (NetCDF) files will reside in (if not given a suitable default is applied)."
    echo
    printf "  %-20s%-s\n" "--dbg" \
	   "whether to compile code(s) with/without debug flags (dbg=1/dbg=0)."
    echo
    echo "------------------------------"
} #usage


##################################################
#
#--        c o m m a n d   l i n e
#
while [ $# -gt 0 ]
do
    arg="$1"
    case $arg in
      --vanilla)
	is_vanilla=1
	;;
      --org)
	tst_org=1
	;;
      --runsimobs)
	tst_runsimobs=1
	;;
      --twin)
	tst_twin=1
	;;
      --jacfwv)
	tst_jacfwv=1
	;;
      --assim_obs)
	shift
	assim_obs=$1
	;;
      --nrlib)
	shift
	external_nrlib=$1
	;;
      --initial_yyyymmdd)
	shift
	initial_yyyymmdd=$1
	;;
      --ndays)
	shift
	ndays=$1
	;;
      --nlp)
	shift
	nlp=$1
	;;
      --filename_prefixes)
	shift
	filename_prefixes="${filename_prefixes} $1"
	shift
	filename_prefixes="${filename_prefixes} $1"
	shift
	filename_prefixes="${filename_prefixes} $1"
	;;
      --pert)
	shift
	pert=$1
	;;
      --runc)
	shift
	runc=$1
	;;
      --disable_prior)
	disable_prior=1
	;;
      --expdir)
        shift
        exp_dir=$1
        ;;
      --outdir)
        shift
        out_dir=$1
        ;;
      --dbg)
	shift
	dbg=$1
	;;
      -h)
        usage
        exit 0
        ;;
      --help)
        usage
        exit 0
        ;;
      -*)
        echo "  ERROR::ignore unexpected option ->>>${arg}<<<-"
        ;;
      *)
        echo "  ERROR::unexpected argument ->>>${arg}<<<-"
        exit 1
        ;;
    esac
    shift
done


##################################################
#
#--        s t a r t   s c r i p t
#
if [ "x${exp_dir}" = "x" ]
then
  cur_dir=`pwd`
  mkdir -p testarea && exp_dir=`mktemp -d -p ${cur_dir}/testarea`
  if [ $? -ne 0 ]
  then
    echo "${script_name}::FATAL::exp_dir could not be generated!"
    exit 1
  fi
fi
if [ "x${out_dir}" = "x" ]
then
  out_dir=${exp_dir}/output
fi
src_dir=${exp_dir}/src
#-- logging
echo "${script_name}::INFO:basedir     ***${basedir}***"
echo "${script_name}::INFO:exp_dir     ***${exp_dir}***"
echo "${script_name}::INFO:src_dir     ***${src_dir}***"
echo "${script_name}::INFO:out_dir     ***${out_dir}***"

#-- create experiment/output directory

mkdir -p ${src_dir}
mkdir -p ${out_dir}

#-- prefer fully qualified path (eases switching between directories)
exp_dir=$(readlink -f ${exp_dir})
out_dir=$(readlink -f ${out_dir})
src_dir=$(readlink -f ${src_dir})

#-- copy Makefile
cp -t ${src_dir} ${bld_dir}/Makefile ${bld_dir}/Makefile.org
#-- module dependency script
cp ${bld_dir}/Mkdepends ${src_dir} && chmod u+x ${src_dir}/Mkdepends
#-- derivative check namelist
if [ -f ${bld_dir}/chkdv.par ]
then
  cp ${bld_dir}/chkdv.par ${exp_dir}
fi


#-- copy files from source directory to compile (i.e. obj directory)
src_files=$(find ${model_dir} -maxdepth 1 -type f)
cp -t ${src_dir}  ${src_files}
esmfwrf_files=$(find ${model_dir}/esmf_wrf_timemgr -maxdepth 1 -type f)
cp -t ${src_dir} ${esmfwrf_files}
#-- iLab environment including derivatives
ilab_files=$(find ${model_dir}/iLab -maxdepth 1 -type f)
cp -t ${src_dir} ${ilab_files}
# AD support files (push/pop mechanism)
cp -r ${model_dir}/iLab/adsupportlib ${src_dir}
# Numerical Recipes library
cp -r ${model_dir}/iLab/nr ${src_dir}
if [ "x${external_nrlib}" != "x" ]
then
  cp ${external_nrlib} ${src_dir}
fi
#-- iLab driver
ilab_driver=$(find ${model_dir}/iLab/driver -maxdepth 1 -type f)
cp -t ${src_dir} ${ilab_driver}



#===================
#
#    build namelist for BEPS run
#
echo "${script_name}::INFO: building BEPS namelist..."
echo "${script_name}::  yyyymmdd_start=${yyyymmdd_start} ndays=${ndays} nlp=${nlp}"
echo "${script_name}::  filename_prefixes ***${filename_prefixes}***"
cd ${exp_dir}
#-- NOTE::output directory must either be
#         - an absolute path (maybe we must be careful with its length...)
#         - a qualified path relative to exp_dir
nml_builder=${basedir}/bld/beps_build_nml.sh
if [ "x${filename_prefixes}" != "x" ]
then
  ${nml_builder} -i_yyyymmdd ${initial_yyyymmdd} -ndays ${ndays} -nlp ${nlp} -filename_prefixes ${filename_prefixes} -iptdir $ipt_dir -outdir $out_dir
else
  ${nml_builder} -i_yyyymmdd ${initial_yyyymmdd} -ndays ${ndays} -nlp ${nlp} -iptdir $ipt_dir -outdir $out_dir
fi
echo "${script_name}::INFO: ...namelist DONE."


#===================
#
#    vanilla
#
if [ ${is_vanilla} -eq 1 ]
then
  echo "${script_name}::INFO:new experiment directory ***${exp_dir}*** is prepared."
  echo "${script_name}::     switch to ${src_dir} in order to compile a selected driver."
  echo "${script_name}::     a compiled driver has to be invoked from within ${exp_dir}."
  exit 0
fi


#===================
#
#    tst_org
#
if [ ${tst_org} = 1 ]
then
  cd ${src_dir}
  #-- create dependency file
  make dep
  #--
  prg=runbepsorg
  exe=${prg}.x
  echo "${script_name}::INFO: compile executable +++${exe}+++ for verification test..."
  make ${exe} DBG=${dbg} 2>&1 | tee ${exp_dir}/compile_${prg}.log
  if [ -f ${exe} ]
  then
    echo "${script_name}::INFO: ...compilation DONE."
  else
    echo "${script_name}::FATAL: compilation failed, exectuable ***${exe}*** was *not* generated"
    exit 1
  fi
  cd ${exp_dir}
  echo "${script_name}::INFO: start BEPS verification run..."
  ${src_dir}/${exe} 2>&1 | tee run_${prg}.log
  rc=$?
  if [ ${rc} -ne 0 ]
  then
    echo "${script_name}::FATAL:BEPS run terminated with non-zero exit code"
  else
    echo "${script_name}::...BEPS run terminated!"
  fi
  exit ${rc}
fi

#===================
#
#    tst_runsimobs
#
if [ ${tst_runsimobs} -eq 1 ]
then
  cd ${src_dir}
  #-- create dependency file
  make dep
  #--
  prg=runsimobs
  exe=${prg}.x
  make ${exe} DBG=${dbg} 2>&1 |& tee ${exp_dir}/compile_${prg}.log
  if [ -f ${exe} ]
  then
    echo "${script_name}::INFO: ...compilation DONE."
  else
    echo "${script_name}::FATAL: compilation failed, exectuable ***${exe}*** was *not* generated"
    exit 1
  fi
  cd ${exp_dir}
  echo "${script_name}::INFO: start BEPS verification run..."
  logfile=run_${prg}.log
  ${src_dir}/${exe} --runc ${runc} 2>&1 | tee ${logfile}
  rc=$?; \
  if [ ${rc} -ne 0 ]
  then
    echo "${script_name}::FATAL:BEPS run terminated with non-zero exit code"
  else
    echo "${script_name}::...BEPS run terminated!"
  fi
  exit ${rc}
fi


#===================
#
#    tst_twin
#
if [ ${tst_twin} -eq 1 ]
then
  cd ${src_dir}
  #-- create dependency file
  make dep
  #-- create simulation driver
  prg=runsimobs
  exe=${prg}.x
  make ${exe} DBG=${dbg} 2>&1 |& tee ${exp_dir}/compile_${prg}.log
  if [ -f ${exe} ]
  then
    echo "${script_name}::INFO: ...compilation DONE."
  else
    echo "${script_name}::FATAL: compilation failed, exectuable ***${exe}*** was *not* generated"
    exit 1
  fi
  #-- create inversion driver
  prg=rundfprv
  exe=${prg}.x
  make ${exe} DBG=${dbg} 2>&1 |& tee ${exp_dir}/compile_${prg}.log
  if [ -f ${exe} ]
  then
    echo "${script_name}::INFO: ...compilation DONE."
  else
    echo "${script_name}::FATAL: compilation failed, exectuable ***${exe}*** was *not* generated"
    exit 1
  fi
  #-- switch to experiment directory
  cd ${exp_dir}
  #-- create pseudo observations
  echo "${script_name}::INFO: start generating pseudo-observations..."
  prg=runsimobs
  exe=${src_dir}/${prg}.x
  logfile=run_${prg}.log
  ${exe} --runc ${runc} 2>&1 | tee ${logfile}
  if [ ! -f obs.nc ]
  then
    echo "${script_name}::FATAL:pseudo observations 'obs.nc' have not been generated"
    exit 1
  else
    echo "${script_name}::...BEPS run terminated and 'obs.nc' was generated!"
  fi
  #-- run inversion
  drvargs="--pert ${pert}"
  if [ ${disable_prior} -eq 1 ]
  then
    drvargs="${drvargs} --disable_prior"
  fi
  echo "${script_name}::INFO: start running inversion with arguments *****${drvargs}*****..."
  prg=rundfprv
  exe=${src_dir}/${prg}.x
  logfile=run_${prg}.log
  ${exe} ${drvargs} 2>&1 | tee ${logfile}
  echo "${script_name}::...inversion has finished (logfile=${logfile})"
fi


#===================
#
#    tst_jacfwv
#
if [ ${tst_jacfwv} = 1 ]
then
  cd ${src_dir}
  #-- create dependency file
  make dep
  #-- create derivative driver
  prg=runderivfwv
  exe=${prg}.x
  make ${exe} DBG=${dbg} 2>&1 |& tee ${exp_dir}/compile_${prg}.log
  if [ -f ${exe} ]
  then
    echo "${script_name}::INFO: ...compilation DONE."
  else
    echo "${script_name}::FATAL: compilation failed, exectuable was *not* generated"
    exit 1
  fi
  cd ${exp_dir}
  echo "${script_name}::INFO: start BEPS vector-forward-derivative Jacobian creation..."
  ${src_dir}/${exe} --mode jacfwv 2>&1 | tee run_${prg}.jacfwv.log
  rc=$?
  echo "${script_name}::INFO:...vector-forward-derivative Jacbobian DONE. (rc=${rc})"
fi


#===================
#
#    assim_obs
#
if [ "x${assim_obs}" != "x" ]
then
  if [ ! -f ${assim_obs} ]
  then
    echo "${script_name}:: observational data file ***${assim_obs}*** not found on system."
    exit 1
  fi
  cd ${src_dir}
  #-- create dependency file
  make dep
  #-- create inversion driver
  prg=rundfprv
  exe=${prg}.x
  make ${exe} DBG=${dbg} 2>&1 |& tee ${exp_dir}/compile_${prg}.log
  if [ -f ${exe} ]
  then
    echo "${script_name}::INFO: ...compilation DONE."
  else
    echo "${script_name}::FATAL: compilation failed, exectuable ***${exe}*** was *not* generated"
    exit 1
  fi
  #-- switch to experiment directory
  cd ${exp_dir}
  #-- link observations
  echo "linking observational file ${assim_obs} to obs.nc in experiment directory"
  ln -s ${assim_obs} obs.nc
  #-- run inversion
  echo "${script_name}::INFO: start running inversion now..."
  prg=rundfprv
  exe=${src_dir}/${prg}.x
  logfile=run_${prg}.log
  ${exe} --pert ${pert} 2>&1 | tee ${logfile}
  echo "${script_name}::...inversion has finished (logfile=${logfile})"
fi
