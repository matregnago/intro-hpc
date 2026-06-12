#
#  @file fix_doxygen_date.sh
#
#  @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @date 2024-02-18
#
# This script fix the date doxygen markup.
#
#!/bin/sh
header=1

print_header()
{
    if [ $header -ne 0 ]
    then
        echo "------ $1 --------"
        header=0
    fi
}

#
# Update the @date with the creation date (first commit date)
#
fix_header_date()
{
    header=1
    date_in_file=`cat $1 | grep -o "@date [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]" | cut  -c7-17`
    date_fix=`git log --date=format:'%Y-%m-%d' -- $1 | grep Date | sort | cut -c9-19 | head -n1`
    if [ -z "$date_in_file" ]
    then
        echo "$1 date is missing and should be $date_fix"
        sed -i -e "s#\(@author.*\)#\1\n * @date $date_fix#" $1
    else
        echo "$1 $date_in_file $date_fix"
        if [ "$date_in_file" != "$date_fix" ]
        then
            sed -i -e "s#@date [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]#@date $date_fix#" $1
        fi
    fi
}

#
# Loop over files
#
files=`git ls-files | grep -v "^\." | grep -v ".*\.md" | grep -v LICENSE | grep -v ".*\.cmake" | grep -v "testing/lin" | grep -v doc/ | grep -v CTest | grep -v cblas.h | grep -v lapacke.h | grep -v "simucore/perfmodels/\.starpu" | grep -v "\.org"`

if [ $# -gt 0 ]
then
    files=$*
fi

for f in $files
do
    #echo $f
    if [ -d $f ]
    then
        continue;
    fi

    fix_header_date $f
done
