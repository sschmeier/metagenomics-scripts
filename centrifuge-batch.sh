#!/usr/bin/env bash

# ##################################################
#
version="0.0.1"              # Sets version variable
#
# HISTORY:
#
# * DATE - v0.0.1  - First Creation
#
# ##################################################

# DO things HERE
function mainScript() {
    echo -n
    #info ${args[@]}
    
    # check if right number of arguments
    if test ${#args[@]} -ne 3; then
        usage
        die "Script expects 3 arguments."
    fi

    IDX=${args[0]}
    INDIR=${args[1]}
    OUTDIR=${args[2]}

    # fasta files instead of fastq?
    optfa=""
    if $fasta ; then
        files="*.fa"
        optfa="-f"
    fi

    if $gzipped ; then
        files=$files".gz"
    fi
    
    filelist=(${INDIR}/$files)
    NUMBER=${#filelist[@]}
    if $combine ; then
        filelist=$(printf ",%s" "${filelist[@]}")
        filelist=${filelist:1}
        #info $filelist
        header "PROCESSING $NUMBER FILES:" $filelist
        info "RUN CENTRIFUGE..."
        centrifuge --quiet -p $processes $optfa -x $IDX -U $filelist -S $OUTDIR/results.txt --report-file $OUTDIR/report.txt
        info "CONVERT REPORT..."
        centrifuge-kreport --show-zeros -x $IDX $OUTDIR/results.txt > $OUTDIR/report-krakenlike.txt
   
        if $krona; then
            info "RUN KRONA..."
            tmpfile=${tmpDir}-krona
            cat $OUTDIR/results.txt | cut -f 1,3 > ${tmpfile}
            ktImportTaxonomy ${tmpfile} -o $OUTDIR/taxonomy.krona.html
        fi
        info 'Compress results'
        gzip $OUTDIR/results.txt
        gzip $OUTDIR/report-krakenlike.txt
        gzip $OUTDIR/report.txt
        info "DONE"
        
    else
        COUNTER=1
        for file in ${INDIR}/$files; do
            header "PROCESSING FILE $COUNTER/$NUMBER:" $file
            OUT=${OUTDIR}/$(basename $file)
            mkdir -p $OUT;
            info "RUN CENTRIFUGE..."
            centrifuge --quiet -p $processes $optfa -x $IDX -U $file -S $OUT/results.txt --report-file $OUT/report.txt
            info "CONVERT REPORT..."
            centrifuge-kreport --show-zeros -x $IDX $OUT/results.txt > $OUT/report-krakenlike.txt
            
            if $krona; then
                info "RUN KRONA..."
                tmpfile=${tmpDir}-krona
                cat $OUT/results.txt | cut -f 1,3 > ${tmpfile}
                ktImportTaxonomy ${tmpfile} -o $OUT/taxonomy.krona.html
            fi
            let COUNTER=COUNTER+1
            
            info 'Compress results'
            gzip $OUTDIR/results.txt
            gzip $OUTDIR/report-krakenlike.txt
            gzip $OUTDIR/report.txt
        done
        info "DONE"
    fi
}

function trapCleanup() {
  echo ""
  # Delete temp files, if any
  if [ -d "${tmpDir}" ] ; then
    rm -r "${tmpDir}"
  fi
  die "Exit trapped. In function: '${FUNCNAME[*]}'"
}

function safeExit() {
  # Delete temp files, if any
  if [ -d "${tmpDir}" ] ; then
    rm -r "${tmpDir}"
  fi
  trap - INT TERM EXIT
  exit
}

# Logging & Feedback
# -----------------------------------------------------
function _alert() {
  if [ "${1}" = "error" ]; then local color="${bold}${red}"; fi
  if [ "${1}" = "warning" ]; then local color="${red}"; fi
  if [ "${1}" = "success" ]; then local color="${green}"; fi
  if [ "${1}" = "debug" ]; then local color="${purple}"; fi
  if [ "${1}" = "header" ]; then local color="${bold}${tan}"; fi
  if [ "${1}" = "input" ]; then local color="${bold}"; fi
  if [ "${1}" = "info" ] || [ "${1}" = "notice" ]; then local color=""; fi
  # Don't use colors on pipes or non-recognized terminals
  if [[ "${TERM}" != "xterm"* ]] || [ -t 1 ]; then color=""; reset=""; fi

  # Print to console when script is not 'quiet'
  if ${quiet}; then return; else
   echo -e "$(date +"%H:%M:%S") ${color}$(printf "[%7s]" "${1}") ${_message}${reset}";
  fi

  # Print to Logfile
  
  if ${printLog} && [ "${1}" != "input" ]; then
    color=""; reset="" # Don't use colors in logs
    echo -e "$(date +"%Y%m%d-%H:%M:%S") $(printf "[%7s]" "${1}") ${_message}" >> "${logFile}";
  fi
}

function die ()       { local _message="${*} Exiting."; echo -e "$(_alert error)"; safeExit;}
function error ()     { local _message="${*}"; echo -e "$(_alert error)"; }
function warning ()   { local _message="${*}"; echo -e "$(_alert warning)"; }
function notice ()    { local _message="${*}"; echo -e "$(_alert notice)"; }
function info ()      { local _message="${*}"; echo -e "$(_alert info)"; }
function debug ()     { local _message="${*}"; echo -e "$(_alert debug)"; }
function success ()   { local _message="${*}"; echo -e "$(_alert success)"; }
function input()      { local _message="${*}"; echo -n "$(_alert input)"; }
function header()     { local _message="== ${*} ==  "; echo -e "$(_alert header)"; }
function verbose()    { if ${verbose}; then debug "$@"; fi }



# Set Base Variables
# ----------------------
scriptName=$(basename "$0")

# Set default Flags
quiet=false
printLog=false
verbose=false
strict=false
debug=false
fasta=false
files="*.fastq"
processes=4
combine=false
gzipped=false
krona=false
args=()

# Set Colors
bold=$(tput bold)
reset=$(tput sgr0)
purple=$(tput setaf 171)
red=$(tput setaf 1)
green=$(tput setaf 76)
tan=$(tput setaf 3)
blue=$(tput setaf 38)
underline=$(tput sgr 0 1)

# Set Temp Directory
tmpDir="/tmp/${scriptName}.$RANDOM.$RANDOM.$RANDOM.$$"
(umask 077 && mkdir "${tmpDir}") || {
  die "Could not create temporary directory! Exiting."
}

# Logging
# -----------------------------------
# Log is only used when the '-l' flag is set.
logFile="${scriptName}-$(date +"%Y%m%d-%H%M%S").log"

# Options and Usage
# -----------------------------------
usage() {
  echo -n "${scriptName} [OPTION] CF-IDX INDIR OUTDIR

Script will look for *.fastq files in INDIR and run centrifuge with the CF-IDX 
index for each file (or combined if -c is set). Results will be put in OUTDIR.

 ${bold}General Options:${reset}
  -q, --quiet       Quiet (no output)
  -l, --log         Print log to file
  -s, --strict      Exit script with null variables.  i.e 'set -o nounset'
  -v, --verbose     Output more information. (Items echoed to 'verbose')
  -d, --debug       Runs script in BASH debug mode (set -x)
  -h, --help        Display this help and exit
      --version     Output version information and exit

 ${bold}Specific Options:${reset}
  -f, --fasta       Query input files are (multi-)FASTA .fa
      --gzip        Input is gzipped compressed.
  -c, --combine     Combine the results of all files 
                    (Much faster, however input-file info lost)
  -k, --krona       Run krona on the centrifuge results
  -p NUM            Number of processes to use [default=4]

"
}

# Iterate over options breaking -ab into -a -b when needed and --foo=bar into
# --foo bar
optstring=h
unset options
while (($#)); do
  case $1 in
    # If option is of type -ab
    -[!-]?*)
      # Loop over each character starting with the second
      for ((i=1; i < ${#1}; i++)); do
        c=${1:i:1}

        # Add current char to options
        options+=("-$c")

        # If option takes a required argument, and it's not the last char make
        # the rest of the string its argument
        if [[ $optstring = *"$c:"* && ${1:i+1} ]]; then
          options+=("${1:i+1}")
          break
        fi
      done
      ;;

    # If option is of type --foo=bar
    --?*=*) options+=("${1%%=*}" "${1#*=}") ;;
    # add --endopts for --
    --) options+=(--endopts) ;;
    # Otherwise, nothing special
    *) options+=("$1") ;;
  esac
  shift
done

set -- "${options[@]}"
unset options

# Print help if no arguments were passed.
# -------------------------------------
[[ $# -eq 0 ]] && set -- "--help"

# Read the options and set stuff
while [[ $1 = -?* ]]; do
  case $1 in
    -h|--help) usage >&2; safeExit ;;
    --version) echo "$(basename $0) ${version}"; safeExit ;;
    -v|--verbose) verbose=true ;;
    -l|--log) printLog=true ;;
    -q|--quiet) quiet=true ;;
    -s|--strict) strict=true;;
    -d|--debug) debug=true;;
    -f|--fasta) fasta=true;;
    -c|--combine) combine=true;;
    -k|--krona) krona=true;;
    --gzip) gzipped=true;;
    -p) shift; processes=$1;;
    --endopts) shift; break ;;
    *) die "invalid option: '$1'." ;;
  esac
  shift
done

# Store the remaining part as arguments.
args+=("$@")

# Trap bad exits with your cleanup function
trap trapCleanup EXIT INT TERM

# Set IFS to preferred implementation
IFS=$' \n\t'

# Exit on error. Append '||true' when you run the script if you expect an error.
set -o errexit

# Run in debug mode, if set
if ${debug}; then set -x ; fi

# Exit on empty variable
if ${strict}; then set -o nounset ; fi

# Bash will remember & return the highest exitcode in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`, for example.
set -o pipefail

# Run your script
mainScript

# Exit cleanly
safeExit

