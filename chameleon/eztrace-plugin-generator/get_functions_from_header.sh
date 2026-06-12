#!/bin/bash
###
#
#  @file get_functions_from_header.sh
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
output_file=$2
header_name=$(echo $header_file |sed 's/\.h//')
if [ ! -f "$header_file" ];then
   echo "'$header_file' is not a file " >&2
   exit 1
fi

tmpdir=$(mktemp -d)
prepoc_header=$tmpdir/preproc.h
processed_header=$tmpdir/processed.h
indent_header_file=$tmpdir/processed_indented.h

function cleanup()
{
    echo "Cleaning up $tmpdir" >&2
    rm -rf "$tmpdir"
}

# run the preprocessor to expand macros
#gcc -E "$header_file"  -o "$prepoc_header"  $CFLAGS
#if [ $? -ne 0 ]; then
#    echo "running 'gcc -E \"$header_file\" -o \"$prepoc_header\" $CFLAGS' failed !" >&2
#    echo -e "You can try to set CFLAGS to fix the problem." >&2
#    exit 1
#fi

cat $header_file > $prepoc_header

# remove all the code that was #include'd
grep -v -P '#\s+\d+\s+"<.+>' $prepoc_header | \
    awk '
{
  if (/^#\s+[0-9]+\s+".+".+1/) {
# when matching a line like this:
# # 1 "/usr/lib/gcc/x86_64-linux-gnu/13/include/stddef.h" 1 3 4
# -> the 1 flag indicates that we entered a file

    ignore_lines++;

  } else if (/^#\s+[0-9]+\s+".+".+2/) {
# when matching a line like this:
# # 17 "netcdf/netcdf.h" 2
# -> the 2 flag indicates that we leaved a file

    ignore_lines--;

  } else if(ignore_lines <= 1) {
# we are in the main file, simply print the line
    print($0);
  }
}'  > "$processed_header"  || exit 1


# make sure everything is on one line
indent -l5000 "$processed_header" -o "$indent_header_file"

# only keep the function prototypes
sed -i '/_z_h_/d' $indent_header_file
sed -i '/_zc_h_/d' $indent_header_file
sed -i '/_ztile_h_/d' $indent_header_file
sed -i '/_zctile_h_/d' $indent_header_file
gcc -fpreprocessed -dD -E -P $indent_header_file -o $2

cleanup
