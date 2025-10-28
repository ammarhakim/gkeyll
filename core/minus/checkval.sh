#!/bin/sh

RED='\033[0;31m'
GREEN='\033[0;32m'

file=$1
out=`cat ${file}_val_err | grep "no leaks are possible"`
out_definitely=`cat ${file}_val_err | grep "definitely lost: 0 bytes in 0 blocks"`
out_indirectly=`cat ${file}_val_err | grep "indirectly lost: 0 bytes in 0 blocks"`
out_possibly=`cat ${file}_val_err | grep "possibly lost: 0 bytes in 0 blocks"`
out_err=`cat ${file}_val_err | grep "ERROR SUMMARY: 0 errors"`

if [ "$out" = "" ]; then # output Does NOT contain: "no leaks are possible"
    if [ "$out_err" = "" ] || [ "$out_definitely" = "" ] || [ "$out_indirectly" = "" ] || [ "$out_possibly" = "" ]; then
        echo "${RED}Unit test" $file "has error issues or memory leaks!"
    else
        echo "${GREEN}Unit test" $file "is valgrind clean"
    fi
else # output contains: "no leaks are possible"
    if [ "$out_err" = "" ]; then # output Does NOT contain: "ERROR SUMMARY: 0 errors"
        echo "${RED}Unit test" $file "has error issues or memory leaks!"
    else # output contains: "ERROR SUMMARY: 0 errors"
        echo "${GREEN}Unit test" $file "is valgrind clean"
    fi
fi
