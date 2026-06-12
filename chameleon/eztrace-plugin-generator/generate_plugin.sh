#!/bin/bash
###
#
#  @file generate_plugin.sh
#  @copyright 2013-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Brieuc Nicolas
#  @date 2025-12-19
#
#  Modified from eztrace-plugin-generator
#  https://gitlab.com/eztrace/eztrace-plugin-generator
###

header_file=$1
if [[ ! -f "$header_file" ]];then
   echo "'$header_file' is not a file " >&2
   exit 1
fi

tmpdir=$(mktemp -d)
prefix=$(realpath $(dirname $0))
template_dir="$prefix/plugin_template"

output_dir=""
PLUGIN_NAME=""
CALLBACK_VARIABLES=""
WRAPPER_FUNCTIONS=""
INTERCEPTIONS=""
PRECISION=""

coreblas_include="#include \"coreblas.h\""
gpucublas_include="#include \"gpucublas.h\""

echo $header_file
if [[ "$header_file" =~ gpucublas ]]; then
    INCLUDE=$gpucublas_include
else
    INCLUDE=$coreblas_include
fi

if [[ "$header_file" =~ \_zc\.h|\_zctile\.h ]]; then
    PRECISION="mixed zc -> zc ds"
elif [[ "$header_file" =~ \_z\.h|\_ztile\.h ]]; then
    PRECISION="normal z -> z c d s"
else
    echo "no match"
fi


function set_plugin_name()
{
    filename=$(basename $header_file)
    PLUGIN_NAME=$(echo $filename |sed 's/\.h//')
    BLAS_NAME=$(echo $PLUGIN_NAME| sed -E 's/_z|_zc|_ztile|_zctile//')

    if [[ "$PLUGIN_NAME" =~ coreblas_zctile ]]; then
        PLUGIN_NAME="coreblas_tcore_zc"
    elif [[ "$PLUGIN_NAME" =~ coreblas_zc ]]; then
        PLUGIN_NAME="coreblas_core_zc"
    elif [[ "$PLUGIN_NAME" =~ coreblas_ztile ]]; then
        PLUGIN_NAME="coreblas_tcore_z"
    elif [[ "$PLUGIN_NAME" =~ coreblas_z ]]; then
        PLUGIN_NAME="coreblas_core_z"
    fi
    output_dir=$(dirname $header_file)/../../eztrace_module
    mkdir -p "$output_dir"

    echo "Creating plugin '$PLUGIN_NAME'"
}

function cleanup()
{
    echo "Cleaning up $tmpdir" >&2
    rm -rf "$tmpdir"
    rm ${PLUGIN_NAME}_clean.h
}


function get_parameters_call()
{
    local parameters="$@,"
    preq=0

    while read  -d',' item; do
	if [ "$preq" = "1" ]; then
	    echo -n ", "
	fi
	if [ "$item" = "..." ]; then
	    return -1
	fi
	if [ "$item" != "void" ]; then
	    preq=1
	    echo -n "$item" |sed 's/.*[ *]//g' |tr -d '[]'
	fi
    done <<< $parameters
}

function process_function()
{
    prototype="$@"
    if [[ "$prototype" =~ ^[[:space:]]*(.+)[[:space:]]+([^ ]+)[[:space:]]*'('(.*)')'[[:space:]]*';'[[:space:]]*$ ]];
    then
	return_type="${BASH_REMATCH[1]}"
	function_name="${BASH_REMATCH[2]}"
	parameters="${BASH_REMATCH[3]}"

	echo "Processing $function_name"
	return_decl="    ";
	return_call="";

	if [ $return_type != "void" ]; then
	    return_decl="    $return_type ret = ";
	    return_call="    return ret;"
	fi

	parameters_call=$(get_parameters_call $parameters)
	if [ $? -ne 0 ]; then
	    echo -e "\tCannot instrument function $function_name"
	    echo -e "\t\tfunction prototype: $prototype"
	    return
	fi

	CALLBACK_VARIABLES="$CALLBACK_VARIABLES$return_type (*lib$function_name)($parameters) = NULL;\n"
	INTERCEPTIONS="${INTERCEPTIONS}INTERCEPT3(\"$function_name\", lib$function_name)\n"

	WRAPPER_FUNCTIONS=$(cat <<EOF | sed -z 's/\n/\\n/g'
$WRAPPER_FUNCTIONS
$return_type $function_name($parameters) {
    FUNCTION_ENTRY;
${return_decl}lib$function_name($parameters_call);
    FUNCTION_EXIT;
$return_call
}
EOF

			 )
    else
        if [[ "$prototype" =~ ^#if[[:space:]]+defined || "$prototype" =~ ^#endif ]]; then
            CALLBACK_VARIABLES="${CALLBACK_VARIABLES}$prototype\n"
            INTERCEPTIONS="${INTERCEPTIONS}$prototype\n"
            WRAPPER_FUNCTIONS="${WRAPPER_FUNCTIONS}$prototype\n"
        else
	        echo "error parsing '$prototype'"
        fi
    fi
}

function process_file()
{
    mkdir -p "$output_dir" 2>/dev/null
    input_file=$1
    output_file=$1
    output_file=$output_dir/$2.c
    echo "Plugin $output_file generated from $input_file"
    sed_script="$tmpdir/script.sed"
    echo 's/@PLUGIN_NAME@/'"$PLUGIN_NAME"'/g' > "$sed_script"
    echo 's/@CALLBACK_VARIABLES@/'"$CALLBACK_VARIABLES"'/g'  >> "$sed_script"
    echo 's/@WRAPPER_FUNCTIONS@/'"$WRAPPER_FUNCTIONS"'/g'  >> "$sed_script"
    echo 's/@INTERCEPTIONS@/'"$INTERCEPTIONS"'/g'  >> "$sed_script"
    echo 's/@INCLUDE@/'"$INCLUDE"'/g'  >> "$sed_script"
    echo 's/@PRECISION@/'"$PRECISION"'/g'  >> "$sed_script"
    echo 's/@NAME@/'"$2"'/g'  >> "$sed_script"

    cat "$input_file" | sed -f "$sed_script"  > "$output_file"
}


set_plugin_name

bash $prefix/get_functions_from_header.sh $header_file ${PLUGIN_NAME}_clean.h
if [ $? -ne 0 ]; then
    echo "get_functions_from_header.sh failed!" >&2
    echo "Try debugging the following command: '$prefix/get_functions_from_header.sh $header_file'"
    cleanup
    exit 1
fi

while read line; do
    process_function "$line"
done < "${PLUGIN_NAME}_clean.h"

for f in "$template_dir/"* ; do
    process_file $f $(echo ${PLUGIN_NAME} |sed 's/coreblas_//')
done

cleanup 2> /dev/null

echo "The plugin was created in $output_dir"
