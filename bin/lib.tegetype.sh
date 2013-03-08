#!/bin/bash
source lib.common.sh
set_explicit_errtrap

NCPU=${NCPU:-4}
if cmd_exists pigz; then
    ZIP=${ZIP:-"pigz -9 -p $NCPU"}
    UNZIP=${UNZIP:-"pigz -d -p $NCPU"}
else
    ZIP=${ZIP:-"gzip -9"}
    UNZIP=${UNZIP:-gunzip}
fi

if cmd_exists pv; then
    PV=${PV:-"pv -f -i 10"}
else
    PV=${PV:-cat}
fi

if [ ! "${BASH_XTRACEFD:-}" ]; then
    XTRACE=${XTRACE:-/dev/null}
    exec {BASH_XTRACEFD}>"$XTRACE"
fi
export BASH_XTRACEFD=$BASH_XTRACEFD

set_ref_var_names () {
    ref_fa=$BASE_DIR/data/ref.$1.fa
    ref_fai=$BASE_DIR/data/ref.$1.fa.fai
}

set_lib_var_names () {
    lib_settings_sh=$BASE_DIR/data/lib.$1.settings.sh
    [ ! -r $lib_settings_sh ] || source $lib_settings_sh
    lib_csv=$BASE_DIR/data/lib.$1.csv
    lib_fa=$BASE_DIR/data/lib.$1.fa
    lib_bt2_idx=$BASE_DIR/data/lib.$ref_name+$1
}

set -x
