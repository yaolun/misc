#!/bin/bash

check_process() {
  echo "$ts: checking $1"
  [ "$1" = "" ]  && return 0
  [ `pgrep -n $1` ] && return 1 || return 0
}  
check_process "chrome"
if [ $? -eq 1 ]; then
    echo "it is running"
fi


# [ $? -eq 1 ] && echo " is running"
# [ $? -eq 0 ] && echo " is not running"
