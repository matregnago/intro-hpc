#
# Check Example lapack_to_chameleon
#

set(TESTLIST
    step0
    step1
    step2
    step3
    step4
    step5
    step6
    step7
    )

foreach(test ${TESTLIST})
    add_test(example_ltm_${test} ./${prec}${test})
endforeach()
